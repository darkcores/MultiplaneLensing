SET(COMMON_SRC
  plummer.cu
  composite.cu
  composite.cpp
  multiplane.cu
  multiplane.cpp
  context.cu
  cpu_context.cpp
  util/cosmology.cu
)

add_library(lens_common SHARED
	${COMMON_SRC}
)

set_target_properties(lens_common PROPERTIES
  CUDA_SEPERABLE_COMPILATION ON)

set_target_properties(lens_common PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(lens_common
  PROPERTIES PUBLIC_HEADER
  "context.h;cpu_context.h")

INSTALL (
  FILES util/cosmology.h util/vector2d.h
  DESTINATION include/multiplanecuda/util
)

# set_property(TARGET lens_common PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
INSTALL(TARGETS lens_common 
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/multiplanecuda
)

# add_executable(lenstest
#  main.cpp
#)

# target_link_libraries(lenstest lens_common)
