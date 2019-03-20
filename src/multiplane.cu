#include "multiplane.h"

#include "util/error.h"
#include <algorithm>
#include <iostream>

Multiplane MultiplaneBuilder::getCuMultiPlane() {
    prepare();

    // Get final lenses from builders
    m_data.clear();
    for (size_t i = 0; i < m_builders.size(); i++) {
        auto lens = m_builders[i].getCuLens();
        PlaneData plane(lens, m_builders[i].redshift());
        m_data.push_back(plane);
    }

    if (m_data.size() == 0 || m_src_data.size() == 0) {
        std::cerr << "No lens and/or source planes given" << std::endl;
        std::terminate();
    }
    cuda = true;
#ifdef __CUDACC__
    // dev_m_data = m_data;
    // dev_m_src_data = m_src_data;
    // plane_ptr = thrust::raw_pointer_cast(&dev_m_data[0]);
    // src_ptr = thrust::raw_pointer_cast(&dev_m_src_data[0]);
    gpuErrchk(cudaMalloc(&plane_ptr, sizeof(PlaneData) * m_data.size()));
    gpuErrchk(cudaMemcpy(plane_ptr, &m_data[0],
                         sizeof(PlaneData) * m_data.size(),
                         cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc(&src_ptr, sizeof(SourcePlane) * m_src_data.size()));
    gpuErrchk(cudaMemcpy(src_ptr, &m_src_data[0],
                         sizeof(SourcePlane) * m_src_data.size(),
                         cudaMemcpyHostToDevice));
#else
    plane_ptr = nullptr;
    src_ptr = nullptr;
#endif
    return Multiplane(m_data.size(), m_src_data.size(), plane_ptr, src_ptr,
                      true);
}

int Multiplane::destroy() {
    if (m_cuda) {
        // So we copy them back to cpu first and then destroy TODO:
        // consider just keeping a vector of them in memory for
        // cleanliness
        size_t psize = m_plane_length * sizeof(PlaneData);
        PlaneData *pptr = (PlaneData *)malloc(psize);
        cpuErrchk(pptr);
        gpuErrchk(cudaMemcpy(pptr, m_plane_ptr, psize, cudaMemcpyDeviceToHost));
        for (int i = 0; i < m_plane_length; i++) {
            pptr[i].lens.destroy();
        }
        free(pptr);
        size_t ssize = m_src_length * sizeof(SourcePlane);
        SourcePlane *sptr = (SourcePlane *)malloc(ssize);
        cpuErrchk(sptr);
        gpuErrchk(
            cudaMemcpy(sptr, m_src_plane_ptr, ssize, cudaMemcpyDeviceToHost));
        for (int i = 0; i < m_src_length; i++) {
            sptr[i].destroy();
        }
        free(sptr);
    } else {
        for (int i = 0; i < m_plane_length; i++) {
            m_plane_data[i].lens.destroy();
        }
        for (int i = 0; i < m_src_length; i++) {
            m_src_data[i].destroy();
        }
    }
    if (m_cuda) {
        gpuErrchk(cudaFree(m_plane_data));
        gpuErrchk(cudaFree(m_src_data));
    } else {
        free(m_plane_data);
        free(m_src_data);
    }
    m_plane_data = NULL;
    m_src_data = NULL;
    return 0;
}

uint8_t Multiplane::traceTheta(Vector2D<float> theta) const {
    // printf("Theta: (%f, %f)\n", theta.x(), theta.y());
    int i_src = 0;
    const uint8_t pixel = 0;
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
    for (int i = 0; i < (m_plane_length - 1); i++) {
        zs = m_plane_ptr[i + 1].redshift;
        // TODO what if sourceplane is before lens
        while (z_src <= zs) {
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
		auto beta = m_plane_ptr[i].lens.getBetaf(theta);
		theta = beta;
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

void Multiplane::traceTheta(Vector2D<float> theta, float *beta_x,
                               float *beta_y, const size_t offset) const {
    int i_src = 0;
    float z_src = m_src_plane_ptr[i_src].redshift();
    // Draw before any lenses first (TODO)
    float zs = z_src + 10000000;
    if (m_plane_length > 0)
        zs = m_plane_ptr[0].redshift;
	// printf("Z init vals: %f and %f\n", z_src, zs);
    while (z_src < zs) {
        // Do source plane(s)
		// printf("Using theta without lens: [ %f ; %f ]\n", theta.x(), theta.y());
		*beta_x = theta.x();
		*beta_y = theta.y();
		beta_x += offset;
		beta_y += offset;
        i_src++;
        if (i_src == m_src_length) {
            // No need to continue now
			return;
        }
        z_src = m_src_plane_ptr[i_src].redshift();
    }

    // Go over each lensplane, and, if we encounter it, source plane.
    for (int i = 0; i < (m_plane_length - 1); i++) {
        zs = m_plane_ptr[i + 1].redshift;
        while (z_src <= zs) {
            // Do source plane(s) that are in front of the current lens
            const float Ds = m_src_plane_ptr[i_src].ds();
            const float Dds = m_src_plane_ptr[i_src].dds();
            const Vector2D<float> s_theta =
                m_plane_ptr[i].lens.getBetaf(theta, Ds, Dds);
			// printf("theta: [ %.8f ; %.8f ]\n", theta.x(), theta.y());
			// printf("s_theta: [ %.8f ; %.8f ]\n", s_theta.x(), s_theta.y());
			*beta_x = s_theta.x();
			*beta_y = s_theta.y();
			beta_x += offset;
			beta_y += offset;
            i_src++;
            if (i_src == m_src_length) {
                // No need to continue now
				return;
            }
            z_src = m_src_plane_ptr[i_src].redshift();
        }
		auto beta = m_plane_ptr[i].lens.getBetaf(theta);
		theta = beta;
    }
	
    // Handle remaining source planes with last lens
    while (i_src < m_src_length) {
        float Ds = m_src_plane_ptr[i_src].ds();
        float Dds = m_src_plane_ptr[i_src].dds();
        Vector2D<float> s_theta =
            m_plane_ptr[m_plane_length - 1].lens.getBetaf(theta, Ds, Dds);
		// printf("final theta: [ %.8f ; %.8f ]\n", theta.x(), theta.y());
		// printf("final s_theta: [ %.8f ; %.8f ]\n", s_theta.x(), s_theta.y());
		*beta_x = s_theta.x();
		*beta_y = s_theta.y();
		beta_x += offset;
		beta_y += offset;
        i_src++;
    }
}