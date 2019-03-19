#include "sources.h"

#include "util/error.h"

#include <iostream>

SourcePlane SourcePlaneBuilder::getCuPlane() {
    cuda = true;
#ifdef __CUDACC__
    // dev_m_points = m_points;
    // ptr = thrust::raw_pointer_cast(&dev_m_points[0]);
    size_t size = sizeof(SourceData) * m_points.size();
	/*
    if (size == 0) {
        std::cerr << "No sources added" << std::endl;
        std::terminate();
    }
	*/
    gpuErrchk(cudaMalloc(&ptr, size));
    gpuErrchk(cudaMemcpy(ptr, &m_points[0], size, cudaMemcpyHostToDevice));
#else
    ptr = nullptr;
#endif
    return SourcePlane(m_redshift, ptr, m_points.size(), true);
}

int SourcePlane::destroy() {
    if (m_points_ptr) {
        if (m_cuda) {
#ifdef __CUDACC__
            gpuErrchk(cudaFree(m_points_ptr));
#endif
        } else {
            free(m_points_ptr);
        }
        m_points_ptr = NULL;
    }
    return 0;
}