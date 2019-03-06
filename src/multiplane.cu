#include "multiplane.h"

#include <algorithm>
#include <iostream>

Multiplane MultiplaneBuilder::getCuMultiPlane() {
    prepare();
#ifdef __CUDACC__
    // dev_m_data = m_data;
    // dev_m_src_data = m_src_data;
    // plane_ptr = thrust::raw_pointer_cast(&dev_m_data[0]);
    // src_ptr = thrust::raw_pointer_cast(&dev_m_src_data[0]);
    cudaMalloc(&plane_ptr, sizeof(PlaneData) * m_data.size());
    cudaMemcpy(plane_ptr, &m_data[0], sizeof(PlaneData) * m_data.size(),
               cudaMemcpyHostToDevice);
    cudaMalloc(&src_ptr, sizeof(SourcePlane) * m_src_data.size());
    cudaMemcpy(src_ptr, &m_src_data[0],
               sizeof(SourcePlane) * m_src_data.size(), cudaMemcpyHostToDevice);
#else
    plane_ptr = nullptr;
    src_ptr = nullptr;
#endif
    return Multiplane(m_data.size(), m_src_data.size(), plane_ptr, src_ptr);
}

void MultiplaneBuilder::cuFree() {
	cudaFree(plane_ptr);
	cudaFree(src_ptr);
}

uint8_t Multiplane::traceTheta(Vector2D<float> theta) const {
    int i_src = 0;
    uint8_t pixel = 0;
    float z_src = m_src_plane_ptr[i_src].redshift();
    // Go over each lensplane, and, if we encounter it, source plane.
#ifndef __CUDA_ARCH__
    std::cout << "Lenses: " << m_plane_length << std::endl;
#endif
    for (int i = 0; i < m_plane_length; i++) {
        float zs = m_plane_ptr[i].redshift;
#ifndef __CUDA_ARCH__
        std::cout << "Theta: " << theta.x() << "," << theta.y() << std::endl;
        std::cout << "Redshifts: " << zs << " src: " << z_src << std::endl;
#endif
        // TODO what if sourceplane is before lens
        while (z_src < zs) {
            // Do source plane(s)
            // TODO check theta with sourceplane
            float Ds = m_src_plane_ptr[i_src].ds();
            float Dds = m_src_plane_ptr[i_src].dds();
            Vector2D<float> s_theta =
                m_plane_ptr[i].lens.getBetaf(theta, Ds, Dds);
            uint8_t p = m_src_plane_ptr[i_src].check_hit(s_theta);
            if (p != 0) {
                // TODO return here or add sources?
                return p;
            }
            i_src++;
            if (i_src == m_src_length) {
                // No need to continue now
                return pixel;
            }
            z_src = m_src_plane_ptr[i_src].redshift();
        }
        if (i < (m_plane_length - 1)) {
#ifndef __CUDA_ARCH__
            std::cout << "Theta updated " << std::endl;
#endif
            auto beta = m_plane_ptr[i].lens.getBetaf(theta);
            theta = beta;
#ifndef __CUDA_ARCH__
            std::cout << "New theta: " << beta.x() << "," << beta.y()
                      << std::endl;
#endif
        }
    }
    // Handle remaining source planes
    while (i_src < m_src_length) {
        float Ds = m_src_plane_ptr[i_src].ds();
        float Dds = m_src_plane_ptr[i_src].dds();
        Vector2D<float> s_theta =
            m_plane_ptr[m_plane_length - 1].lens.getBetaf(theta, Ds, Dds);
        uint8_t p = m_src_plane_ptr[i_src].check_hit(s_theta);
        if (p != 0) {
            // TODO return here or add sources?
            return p;
        }
        i_src++;
    }
    return pixel;
}
