#include "multiplane.h"

#include <algorithm>
#include <iostream>

Multiplane MultiplaneBuilder::getCuMultiPlane() {
    prepare();

    // Get final lenses from builders
    for (size_t i = 0; i < m_builders.size(); i++) {
        auto lens = m_builders[i]->getCuLens();
        PlaneData plane(lens, m_builders[i]->redshift());
        m_data.push_back(plane);
    }

    cuda = true;
#ifdef __CUDACC__
    // dev_m_data = m_data;
    // dev_m_src_data = m_src_data;
    // plane_ptr = thrust::raw_pointer_cast(&dev_m_data[0]);
    // src_ptr = thrust::raw_pointer_cast(&dev_m_src_data[0]);
    cudaMalloc(&plane_ptr, sizeof(PlaneData) * m_data.size());
    cudaMemcpy(plane_ptr, &m_data[0], sizeof(PlaneData) * m_data.size(),
               cudaMemcpyHostToDevice);
    cudaMalloc(&src_ptr, sizeof(SourcePlane) * m_src_data.size());
    cudaMemcpy(src_ptr, &m_src_data[0], sizeof(SourcePlane) * m_src_data.size(),
               cudaMemcpyHostToDevice);
#else
    plane_ptr = nullptr;
    src_ptr = nullptr;
#endif
    return Multiplane(m_data.size(), m_src_data.size(), plane_ptr, src_ptr);
}

void MultiplaneBuilder::cuFree() {
#ifdef __CUDACC__
    cudaFree(plane_ptr);
    cudaFree(src_ptr);
#endif
}

uint8_t Multiplane::traceTheta(Vector2D<float> theta) const {
	// printf("Theta: (%f, %f)\n", theta.x(), theta.y());
    int i_src = 0;
    uint8_t pixel = 0;
    float z_src = m_src_plane_ptr[i_src].redshift();
    // Draw before any lenses first (TODO)
    float zs = z_src + 1;
    if (m_plane_length > 0)
        zs = m_plane_ptr[0].redshift;
    while (z_src < zs) {
        // Do source plane(s)
        uint8_t p = m_src_plane_ptr[i_src].check_hit(theta);
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

    // Go over each lensplane, and, if we encounter it, source plane.
    for (int i = 0; i < m_plane_length; i++) {
        zs = m_plane_ptr[i].redshift;
		/*
		#ifdef __CUDA_ARCH__
		if (threadIdx.x < 3) {
			printf("Theta %d: (%f; %f)\n", threadIdx.x, theta.x(), theta.y());
		}
		#else 
		printf("Theta: (%f; %f)\n", theta.x(), theta.y());
		#endif
		*/
        // TODO what if sourceplane is before lens
        while (z_src < zs) {
            // Do source plane(s)
            // TODO check theta with sourceplane
            float Ds = m_src_plane_ptr[i_src].ds();
            float Dds = m_src_plane_ptr[i_src].dds();
            Vector2D<float> s_theta =
                m_plane_ptr[i - 1].lens.getBetaf(theta, Ds, Dds);
            uint8_t p = m_src_plane_ptr[i_src].check_hit(s_theta);
			/*
#ifdef __CUDA_ARCH__
			if (threadIdx.x < 3) {
				printf("S Theta %d: (%f; %f)\n", threadIdx.x, theta.x(), theta.y());
			}
#else 
			printf("S Theta: (%f; %f)\n", theta.x(), theta.y());
#endif
			*/
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
            auto beta = m_plane_ptr[i].lens.getBetaf(theta);
            theta = beta;
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
