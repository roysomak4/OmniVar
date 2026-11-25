add_executable(omnivar_unit_tests
    main.cpp
    test_cli.cpp
    test_csv_slice_provider.cpp
    test_num_format.cpp
    test_regex_hash.cpp
    test_realign.cpp
    test_reject.cpp
    test_region_builder.cpp
    test_stats_containers.cpp
)
target_include_directories(omnivar_unit_tests SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/tests/third_party)
target_compile_definitions(omnivar_unit_tests PRIVATE OMNIVAR_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(omnivar_unit_tests PRIVATE omnivar_lib)

add_test(NAME unit COMMAND omnivar_unit_tests)
