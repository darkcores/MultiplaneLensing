# cmake_minimum_required (VERSION 3.9.0)
# project (example LANGUAGES CXX C)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
# set(default_build_type "RELEASE")

# include(CheckIncludeFile)
# check_include_file(getopt.h HAVE_GETOPT_H)

# include_directories("../src/")

add_executable(example
  example.cpp
)

add_executable(example_cpu
  example_cpu.cpp
)

target_link_libraries(example
  PRIVATE lens_common
)

target_link_libraries(example_cpu
  PRIVATE lens_common
)
