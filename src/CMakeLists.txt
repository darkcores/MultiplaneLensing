SET(COMMON_SRC
  # lens.cpp
)

add_library(lens_common
	${COMMON_SRC}
)

set_property(TARGET lens_common PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
