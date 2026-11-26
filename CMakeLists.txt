cmake_minimum_required(VERSION 3.25)

# Release builds use -O2, not CMake's default -O3: on the Oncospan panel (2026-11-24, GCC 14, static)
# -O2 was ~3% faster than -O3 on one thread and equal on 16; LTO, -march=x86-64-v3 and PGO gave
# nothing worth their cost. Set before project() so these become
# the cache defaults; -DCMAKE_CXX_FLAGS_RELEASE=... still overrides them.
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "C flags for Release builds")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "C++ flags for Release builds")

project(omnivar
    VERSION 0.1.0
    DESCRIPTION "omnivar: variant caller for targeted and whole-genome sequencing"
    LANGUAGES C CXX)


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(OMNIVAR_STATIC      "Link a fully static binary (Linux) / static htslib + compression libs (macOS)" OFF)
option(OMNIVAR_BUILD_TESTS "Build unit tests" ON)
option(OMNIVAR_WERROR      "Treat compiler warnings as errors" OFF)
option(OMNIVAR_SANITIZE    "Build with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(OMNIVAR_LIBDEFLATE  "Use libdeflate in htslib when available" ON)
set(OMNIVAR_HTSLIB_TARBALL "" CACHE FILEPATH
    "Local htslib release tarball (offline builds); downloaded when empty")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    # Floating-point std::to_chars needs macOS 13.3.
    set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3" CACHE STRING "Minimum macOS version" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compiler checks: GCC >= 14, Clang >= 18, AppleClang >= 16
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
    message(FATAL_ERROR "GCC >= 14 required (found ${CMAKE_CXX_COMPILER_VERSION}); "
                        "pass -DCMAKE_CXX_COMPILER=g++-14")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18)
    message(FATAL_ERROR "Clang >= 18 required (found ${CMAKE_CXX_COMPILER_VERSION})")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16)
    message(FATAL_ERROR "AppleClang >= 16 (Xcode 16) required (found ${CMAKE_CXX_COMPILER_VERSION})")
endif()

# ---------------------------------------------------------------------------
# Warnings, sanitizers, static linking
# ---------------------------------------------------------------------------
add_library(omnivar_options INTERFACE)
target_compile_options(omnivar_options INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
    $<$<BOOL:${OMNIVAR_WERROR}>:-Werror>)

if(OMNIVAR_SANITIZE)
    target_compile_options(omnivar_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined)
    target_link_options(omnivar_options INTERFACE -fsanitize=address,undefined)
    # libc++/libc++abi false positive on every caught exception; see the file.
    target_sources(omnivar_options INTERFACE ${PROJECT_SOURCE_DIR}/src/util/asan_default_options.cpp)
endif()

if(OMNIVAR_STATIC)
    if(OMNIVAR_SANITIZE)
        message(FATAL_ERROR "OMNIVAR_STATIC and OMNIVAR_SANITIZE cannot be combined")
    endif()
    if(APPLE)
        # libSystem stays dynamic; prefer static archives for everything else.
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a .tbd .dylib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
        target_link_options(omnivar_options INTERFACE -static)
    endif()
endif()

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(Threads REQUIRED)
include(cmake/htslib.cmake)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
configure_file(src/version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/generated/version.hpp @ONLY)

add_library(omnivar_lib STATIC
    src/annotate/var2vcf.cpp
    src/annotate/var2vcf_paired.cpp
    src/app/launcher.cpp
    src/core/cigar.cpp
    src/core/cigar_modifier.cpp
    src/core/cigar_parser.cpp
    src/core/realigner.cpp
    src/core/sv.cpp
    src/core/scope.cpp
    src/core/variant.cpp
    src/core/variation_utils.cpp
    src/core/vars_builder.cpp
    src/cli/cli.cpp
    src/io/bam.cpp
    src/io/csv_slice_provider.cpp
    src/io/faidx_provider.cpp
    src/io/record_preprocessor.cpp
    src/io/reference.cpp
    src/modes/pipeline.cpp
    src/modes/simple.cpp
    src/modes/somatic.cpp
    src/region/region_builder.cpp
    src/stats/fisher.cpp
    src/util/num_format.cpp
    src/util/hash_order.cpp
    src/util/regex_engine.cpp
    src/util/number_parse.cpp
)
target_include_directories(omnivar_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(omnivar_lib PUBLIC omnivar_options htslib::htslib)

add_executable(omnivar src/main.cpp)
target_link_libraries(omnivar PRIVATE omnivar_lib)

install(TARGETS omnivar RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(OMNIVAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/unit)

    add_test(NAME cli.version COMMAND omnivar --version)
    set_tests_properties(cli.version PROPERTIES
        PASS_REGULAR_EXPRESSION "^omnivar ${PROJECT_VERSION}, htslib [0-9]")
    add_test(NAME cli.help COMMAND omnivar -H)
    set_tests_properties(cli.help PROPERTIES PASS_REGULAR_EXPRESSION "usage: omnivar")
    add_test(NAME cli.unknown_option COMMAND omnivar --no-such-option)
    set_tests_properties(cli.unknown_option PROPERTIES WILL_FAIL TRUE)

    # Golden parity gates (no JVM): tests/testkit (omnivar-testkit). Extend as milestones land.
    add_subdirectory(tests/testkit)
    if(EXISTS ${PROJECT_SOURCE_DIR}/docs/golden/cases)
        set(TK $<TARGET_FILE:omnivar-testkit>)
        set(BIN --binary $<TARGET_FILE:omnivar>)
        # TSV at every stage: k0_U (no realignment, no SV), U (realignment, M5), full (SV, M6);
        # Simple and Somatic (tumour/normal, M7).
        foreach(stage k0_U U full)
            foreach(mode Simple Somatic)
                string(TOLOWER ${mode} m)
                add_test(NAME golden.${stage}.${m} COMMAND ${TK} golden ${BIN} --stage ${stage} --mode ${mode} -j 4)
                # VCF vs the jar | teststrandbias.R | var2vcf_valid.pl -A chain (paired:
                # testsomatic.R | var2vcf_paired.pl), checked-in expectations; bcftools parse check.
                add_test(NAME golden.vcf.${stage}.${m} COMMAND ${TK} vcf ${BIN} --stage ${stage} --mode ${mode} -j 4)
            endforeach()
        endforeach()
        # -p (M10): adding -p never changes the calls (TSV and VCF, single sample and tumour/normal).
        add_test(NAME golden.reject_invariant.full COMMAND ${TK} reject-invariant ${BIN} -j 4)
        add_test(NAME golden.reject_invariant.U COMMAND ${TK} reject-invariant ${BIN} --extra-args=-U -j 4)
        add_test(NAME golden.reject_invariant.k0_U COMMAND ${TK} reject-invariant ${BIN} "--extra-args=-k 0 -U" -j 4)
        # Worker pool (M9): -th 4 must write exactly what -th 1 writes (golden regions split into chunks).
        add_test(NAME golden.threads COMMAND ${TK} threads ${BIN} -th 4 -j 2)
        # Synthetic rows (many contigs, SVs, clusters) against the R + Perl chain run live.
        find_program(RSCRIPT_EXECUTABLE Rscript)
        find_program(PERL_EXECUTABLE perl)
        if(RSCRIPT_EXECUTABLE AND PERL_EXECUTABLE)
            add_test(NAME golden.vcf.synthetic COMMAND ${TK} vcf ${BIN} --synthetic -k NO_CASES -j 4)
        endif()
    endif()
endif()
