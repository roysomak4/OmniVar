# Pinned htslib, built from the release tarball with the network/plugin features disabled.
# Produces the imported target htslib::htslib (static libhts.a plus its compression libs).

include(ExternalProject)

set(HTSLIB_VERSION 1.22.1)
set(HTSLIB_SHA256 3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7)

# Source tarball, in order: OMNIVAR_HTSLIB_TARBALL (explicit file), then the shared download
# cache, which is filled once from the release URL. The cache is shared by every build directory
# (and by the Docker build through a BuildKit cache mount), so htslib is downloaded once per
# machine. The SHA-256 is checked when the file enters the cache and again by ExternalProject.
set(_hts_name "htslib-${HTSLIB_VERSION}.tar.bz2")
set(_hts_release_url "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/${_hts_name}")

if(DEFINED ENV{OMNIVAR_DOWNLOAD_CACHE})
    set(_cache_default "$ENV{OMNIVAR_DOWNLOAD_CACHE}")
elseif(DEFINED ENV{XDG_CACHE_HOME})
    set(_cache_default "$ENV{XDG_CACHE_HOME}/omnivar/downloads")
elseif(DEFINED ENV{HOME})
    set(_cache_default "$ENV{HOME}/.cache/omnivar/downloads")
else()
    set(_cache_default "${CMAKE_BINARY_DIR}/_downloads")
endif()
set(OMNIVAR_DOWNLOAD_CACHE "${_cache_default}" CACHE PATH
    "Directory that keeps downloaded third-party sources (htslib) across build directories")

if(OMNIVAR_HTSLIB_TARBALL)
    set(_hts_url "${OMNIVAR_HTSLIB_TARBALL}")
else()
    set(_hts_cached "${OMNIVAR_DOWNLOAD_CACHE}/${_hts_name}")
    if(EXISTS "${_hts_cached}")
        file(SHA256 "${_hts_cached}" _hts_have)
        if(NOT _hts_have STREQUAL HTSLIB_SHA256)
            message(WARNING "htslib: ${_hts_cached} has the wrong SHA-256; downloading again")
            file(REMOVE "${_hts_cached}")
        endif()
    endif()
    if(NOT EXISTS "${_hts_cached}")
        message(STATUS "htslib: downloading ${_hts_release_url} into ${OMNIVAR_DOWNLOAD_CACHE}")
        file(MAKE_DIRECTORY "${OMNIVAR_DOWNLOAD_CACHE}")
        # Download to a per-process name, then rename, so parallel configures cannot see a partial file.
        string(RANDOM LENGTH 8 _hts_tmp_suffix)
        set(_hts_tmp "${_hts_cached}.${_hts_tmp_suffix}.part")
        file(DOWNLOAD "${_hts_release_url}" "${_hts_tmp}"
             EXPECTED_HASH SHA256=${HTSLIB_SHA256} STATUS _hts_dl TLS_VERIFY ON)
        list(GET _hts_dl 0 _hts_dl_code)
        if(NOT _hts_dl_code EQUAL 0)
            file(REMOVE "${_hts_tmp}")
            list(GET _hts_dl 1 _hts_dl_msg)
            message(FATAL_ERROR "htslib: download failed (${_hts_dl_msg}). Offline: pass "
                                "-DOMNIVAR_HTSLIB_TARBALL=/path/to/${_hts_name}")
        endif()
        file(RENAME "${_hts_tmp}" "${_hts_cached}")
    else()
        message(STATUS "htslib: using cached ${_hts_cached}")
    endif()
    set(_hts_url "${_hts_cached}")
endif()

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

set(_hts_deflate_flag --without-libdeflate)
set(_hts_deflate_lib "")
if(OMNIVAR_LIBDEFLATE)
    find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h)
    find_library(LIBDEFLATE_LIBRARY NAMES deflate libdeflate)
    if(LIBDEFLATE_INCLUDE_DIR AND LIBDEFLATE_LIBRARY)
        set(_hts_deflate_flag --with-libdeflate)
        set(_hts_deflate_lib ${LIBDEFLATE_LIBRARY})
        message(STATUS "htslib: using libdeflate ${LIBDEFLATE_LIBRARY}")
    else()
        message(STATUS "htslib: libdeflate not found, using zlib only")
    endif()
endif()

# Flags for htslib's own C compile: optimised, and sanitised when the C++ code is.
set(_hts_cflags "-O2 -g")
set(_hts_ldflags "")
if(OMNIVAR_SANITIZE)
    string(APPEND _hts_cflags " -fsanitize=address,undefined -fno-omit-frame-pointer")
    # configure links test programs with CFLAGS objects; they need the sanitizer runtime too.
    string(APPEND _hts_ldflags " -fsanitize=address,undefined")
endif()
if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET)
    string(APPEND _hts_cflags " -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    foreach(_arch IN LISTS CMAKE_OSX_ARCHITECTURES)
        string(APPEND _hts_cflags " -arch ${_arch}")
    endforeach()
endif()

# Include/library search paths for the compression libraries (e.g. Homebrew on macOS).
set(_hts_cppflags "")
foreach(_dir IN LISTS ZLIB_INCLUDE_DIRS BZIP2_INCLUDE_DIR LIBLZMA_INCLUDE_DIRS LIBDEFLATE_INCLUDE_DIR)
    if(_dir)
        string(APPEND _hts_cppflags " -I${_dir}")
    endif()
endforeach()
foreach(_lib IN LISTS ZLIB_LIBRARIES BZIP2_LIBRARIES LIBLZMA_LIBRARIES _hts_deflate_lib)
    if(_lib AND EXISTS "${_lib}")
        get_filename_component(_libdir "${_lib}" DIRECTORY)
        string(APPEND _hts_ldflags " -L${_libdir}")
    endif()
endforeach()

set(_hts_prefix ${CMAKE_BINARY_DIR}/_deps/htslib)
set(_hts_src ${_hts_prefix}/src/htslib_ext)
set(_hts_lib ${_hts_src}/libhts.a)

ExternalProject_Add(htslib_ext
    PREFIX            ${_hts_prefix}
    URL               ${_hts_url}
    URL_HASH          SHA256=${HTSLIB_SHA256}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_IN_SOURCE   TRUE
    CONFIGURE_COMMAND
        ./configure
            CC=${CMAKE_C_COMPILER}
            "CFLAGS=${_hts_cflags}"
            "CPPFLAGS=${_hts_cppflags}"
            "LDFLAGS=${_hts_ldflags}"
            --disable-libcurl --disable-gcs --disable-s3 --disable-plugins
            ${_hts_deflate_flag}
    BUILD_COMMAND     make -j lib-static
    INSTALL_COMMAND   ""
    BUILD_BYPRODUCTS  ${_hts_lib}
    LOG_CONFIGURE     TRUE
    LOG_BUILD         TRUE
    LOG_OUTPUT_ON_FAILURE TRUE)

# The include directory must exist at configure time for the imported target.
file(MAKE_DIRECTORY ${_hts_src})

add_library(htslib::htslib STATIC IMPORTED GLOBAL)
set_target_properties(htslib::htslib PROPERTIES
    IMPORTED_LOCATION ${_hts_lib}
    INTERFACE_INCLUDE_DIRECTORIES ${_hts_src})
target_link_libraries(htslib::htslib INTERFACE
    ${_hts_deflate_lib} ZLIB::ZLIB BZip2::BZip2 LibLZMA::LibLZMA Threads::Threads m)
add_dependencies(htslib::htslib htslib_ext)
