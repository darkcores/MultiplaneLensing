include_directories("../src/")

add_executable(vec2d
  vector2d_tests.cpp
  vector2d_tests.cu
  )

target_link_libraries(vec2d
  ${GTEST_BOTH_LIBRARIES}
)

add_test(AllTestsInVec2D vec2d)
