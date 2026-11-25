# omnivar-testkit: test and validation harness (golden, VCF, threads and -p gates; oracle and
# panel parity against the reference jar; truth-set benchmark). Not installed.
add_executable(omnivar-testkit
    bench.cpp
    common.cpp
    golden.cpp
    main.cpp
    oracle.cpp
    panel.cpp
    parity.cpp
    reject.cpp
    release_notes.cpp
    synthetic.cpp
    threads.cpp
    vcf.cpp
    vcfcmp.cpp
)
target_compile_definitions(omnivar-testkit PRIVATE OMNIVAR_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(omnivar-testkit PRIVATE omnivar_options htslib::htslib Threads::Threads)
