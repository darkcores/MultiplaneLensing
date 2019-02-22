SET(COMMON_SRC
  "util/vector2d.cpp"
  "util/vector2d.cu"
)

add_library(lens_common
	${COMMON_SRC}
)

set_property(TARGET lens_common PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

add_executable(lenstest
  "main.cpp"
)

# target_link_libraries(lenstest lens_common)
