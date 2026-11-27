add_executable(senri_unit_tests
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
target_include_directories(senri_unit_tests SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/tests/third_party)
target_compile_definitions(senri_unit_tests PRIVATE SENRI_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(senri_unit_tests PRIVATE senri_lib)

add_test(NAME unit COMMAND senri_unit_tests)
