#pragma once

#include "composite.h"

#include <cuda_runtime_api.h>
#ifdef __CUDACC__
#include <thrust/device_vector.h>
#endif
#include <thrust/host_vector.h>

class PlanaData {
public:
	PlanaData(CompositeLens lens, double z);

	CompositeLens lens;
	double redshift;
};

class MultiplaneBuilder {
private:
	
};

class Multiplane {
  private:

  public:
    Multiplane();
};
