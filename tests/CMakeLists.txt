include_directories("../src/")

add_executable(unit_tests
  util_tests.cpp
  plummer_tests.cpp
  composite_tests.cpp
  sources_tests.cpp
)

add_executable(cu_unit_tests
  util_tests.cu
  plummer_tests.cu
  composite_tests.cu
)

target_link_libraries(unit_tests
  ${GTEST_BOTH_LIBRARIES}
  lens_common
)

target_link_libraries(cu_unit_tests
  ${GTEST_BOTH_LIBRARIES}
  lens_common
)

add_test(NAME AllTests COMMAND unit_tests)
add_test(NAME AllTestsCuda COMMAND cu_unit_tests)
