cmake_minimum_required (VERSION 3.9.0)
project (LensTestCuda LANGUAGES CUDA CXX C)
set(CMAKE_CXX_COMPILER_WORKS 1)

SET(GCC_FLAGS "-O3 -std=c++14 -Wall -Wextra -pedantic")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_FLAGS}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)
check_include_file(getopt.h HAVE_GETOPT_H)

include(CheckLanguage)
check_language(CUDA)
set(CMAKE_CUDA_FLAGS "-arch=sm_50 -gencode=arch=compute_61,code=sm_61 -dlink --maxrregcount=64 -lineinfo")
# set(CMAKE_CUDA_FLAGS "--default-stream per-thread -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_60,code=sm_60 -Xcompiler -fopenmp")

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

subdirs(src)
subdirs(tests)
