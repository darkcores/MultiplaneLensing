SET(COMMON_SRC
  plummer.cpp
  grid.cpp
  multiplane.cpp
)

add_library(lens_common
	${COMMON_SRC}
)

set_property(TARGET lens_common PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

add_executable(lenstest
  main.cpp
)

# target_link_libraries(lenstest lens_common)
