# omnivar-testkit: test and validation harness (regression cases, thread and -p invariants,
# truth-set benchmark, release notes). Not installed.
add_executable(omnivar-testkit
    bench.cpp
    common.cpp
    main.cpp
    regress.cpp
    reject.cpp
    release_notes.cpp
    threads.cpp
)
target_compile_definitions(omnivar-testkit PRIVATE OMNIVAR_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(omnivar-testkit PRIVATE omnivar_options htslib::htslib Threads::Threads)
