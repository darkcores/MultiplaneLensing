SET(COMMON_SRC
  plummer.cu
  composite.cu
  composite.cpp
  # sources.cu
  # sources.cpp
  # multiplane.cu
  # multiplane.cpp
  # multiplane_cu.cu
  util/cosmology.cu
)

add_library(lens_common
	${COMMON_SRC}
)

# set_property(TARGET lens_common PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# add_executable(lenstest
#  main.cpp
#)

# target_link_libraries(lenstest lens_common)
