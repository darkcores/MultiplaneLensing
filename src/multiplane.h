#pragma once

#include "composite.h"

#include <cuda_runtime_api.h>
#ifdef __CUDACC__
#include <thrust/device_vector.h>
#endif
#include <thrust/host_vector.h>

class PlaneData {
public:
	PlaneData(CompositeLens lens, double z);
	__device__ PlaneData() {}

	CompositeLens lens;
	double redshift;

	bool operator<(const PlaneData &cmp) const {
		return redshift < cmp.redshift;
	}
};

class MultiplaneBuilder {
private:
	thrust::host_vector<PlaneData> m_data;
	
public:
	MultiplaneBuilder();

	void addPlane(CompositeLens &lens, double z);
	
};

class Multiplane {
  private:
	PlaneData *m_plane;

  public:
    Multiplane();

};
