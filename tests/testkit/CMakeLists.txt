# senri-testkit: test and validation harness (regression cases, thread and -p invariants,
# truth-set benchmark, release notes). Not installed.
add_executable(senri-testkit
    bench.cpp
    common.cpp
    main.cpp
    regress.cpp
    reject.cpp
    release_notes.cpp
    threads.cpp
)
target_compile_definitions(senri-testkit PRIVATE SENRI_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(senri-testkit PRIVATE senri_options htslib::htslib Threads::Threads)
