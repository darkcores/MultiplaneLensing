cmake_minimum_required (VERSION 3.9.0)
project (LensTestCuda LANGUAGES CUDA CXX C)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(default_build_type "RELEASE")

SET(GCC_FLAGS "-O3 -std=c++14 -Wall -Wextra -pipe")
if (COVERAGE)
  SET(GCC_FLAGS "-std=c++14 -O0 -Wall -Wextra")
endif(COVERAGE)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_FLAGS}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)
check_include_file(getopt.h HAVE_GETOPT_H)

# include(FindOpenMP)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include(CheckLanguage)
check_language(CUDA)
set(CMAKE_CUDA_FLAGS "-arch=sm_30 -gencode=arch=compute_30,code=sm_30 -dlink --maxrregcount=32")
if (VSC_BUILD) 
  set(CMAKE_CUDA_FLAGS "-arch=sm_60 -gencode=arch=compute_60,code=sm_60 -dlink --maxrregcount=32")
endif(VSC_BUILD)
if (SM_VER)
  set(CMAKE_CUDA_FLAGS "-arch=sm_${SM_VER} -gencode=arch=compute_${SM_VER},code=sm_${SM_VER} -dlink --maxrregcount=32 -lineinfo")
endif(SM_VER)
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
endif(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
# -Xptxas -dlcm=cg
# -ftz=true
# -prec-div=false
# -prec-sqrt=false
# set(CMAKE_CUDA_FLAGS "--default-stream per-thread -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_60,code=sm_60 -Xcompiler -fopenmp")

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories("./src")

if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  find_package(GTest)
  if (GTest_FOUND)
	enable_testing()
	include_directories(${GTEST_INCLUDE_DIRS})
	subdirs(tests)
  endif()
endif()

subdirs(src)
subdirs(example)
