cmake_minimum_required (VERSION 3.9.0)
project (LensTestCuda LANGUAGES CUDA CXX C)
set(CMAKE_CXX_COMPILER_WORKS 1)

SET(GCC_FLAGS "-O3 -std=c++14 -Wall -Wextra -pedantic")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_FLAGS}")

include(CheckIncludeFile)
check_include_file(getopt.h HAVE_GETOPT_H)

include(CheckLanguage)
check_language(CUDA)
# set(CMAKE_CUDA_FLAGS "--default-stream per-thread -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_60,code=sm_60 -Xcompiler -fopenmp")

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

subdirs(src)
subdirs(tests)
