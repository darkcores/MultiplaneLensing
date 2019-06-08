[![pipeline status](http://gitlab.darkcores.net/bapr/lenstest/badges/master/pipeline.svg)](http://gitlab.darkcores.net/bapr/lenstest/commits/master)
[![coverage report](http://gitlab.darkcores.net/bapr/lenstest/badges/master/coverage.svg)](http://gitlab.darkcores.net/bapr/lenstest/commits/master)

# Multiplane Lensing With CUDA

CUDA library for calculating beta vectors for multiplane lensing
situations using multiple Plummer lenses. This is in development and
many problems may occur. This has only been tested on linux (Arch
Linux and CentOS).

# Building

Build requirements:

	g++
	cmake
	cuda
	gtest (optional)

## Default

By default this library will be build for sm `sm_30`. If you have a
newer card, you might get better results specifying a newer arch.

    cmake CMakeLists.txt
    make
	
### GPU Arch optimize

To specify another arch, set `SM_VER` for cmake

	cmake CMakeLists.txt -DSM_VER=sm_61

To build for usage on the [VSC](https://www.vscentrum.be/) (Tesla
p100) use:

	cmake CMakeLists.txt -DVSC_BUILD=1

## Debug

    cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=DEBUG
    make

To build the unit tests you need to have GTest installed. If GTest is
not found, the tests are skipped. The tests are only build in the
DEBUG configuration.

### Running tests (CPU)

    ./bin/unit_tests

### Running tests (CUDA)

    ./bin/cu_unit_tests
	
# Installation

For Arch Linux a PKGBUILD can be found [here](https://gist.github.com/darkcores/d0ee2b9d83c64d2a85aec3773a9ccced)

By default it installs the library to `/usr/local/` on linux, on
archlinux this doesn't work out of the box (`/usr/local/lib` is not in
the ld path, depends on distro).

	sudo make install
	
For archlinux and others you can specify `DESTDIR`.

	sudo make DESTDIR="$pkgdir/" install

Alternatively add to the ld config if you encounter
and error like this:

	./example: error while loading shared libraries: liblens_common.so: cannot open shared object file: No such file or directory
	
By executing the following commands.

	sudo sh -c 'echo "/usr/local/lib" >> /etc/ld.so.conf.d/multiplanecuda.conf'
	sudo ldconfig

# API

To use this library see `example/example.cpp`.

## Documentation

You can run `Doxygen` to generate the documentation. Some
documentation is still missing or outdated.

# Example

For the example there is a notebook that generates some test
parameters. A file generated by this notebook can be used as a test
file.

# References

This project was made for my bachelor's thesis for the School for
Information Technology from the tUL.

[GRALE2](https://github.com/j0r1/GRALE2) has often been used as a
reference for this project.
