#pragma once

#include "plummer.h"
#include "util/error.h"
#include <iostream>
#include <vector>

class CompositeLens {
  private:
    Plummer *__restrict__ m_lenses;
#ifdef __CUDACC__
    float4 *__restrict__ m_lens_int;
#else
    void *m_lens_int;
#endif
    const int m_lenses_size;
    const bool m_cuda;
    float m_mass_sheet;

  public:
    CompositeLens(Plummer *lenses, const int lenses_size, const double Dd = 0.0,
                  const double scale = 1.0, const bool cuda = false);

    int destroy();

    __host__ void update(const float *factors) {
        for (int i = 0; i < m_lenses_size; i++) {
            m_lenses[i].update(factors[i]);
        }
    }

#ifdef __CUDACC__
    __device__ void update(const float &factor, const int idx) {
        m_lenses[idx].update(factor);
        m_lens_int[idx] = m_lenses[idx].f4();
    }

    /**
     * Get alpha vector (single precision, with scaling).
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlpha(const float2 &theta,
                                        const float &mass_scale = 0.0) const {
        const unsigned int DIM = 128;
        float2 alpha{0, 0};

#ifdef __CUDA_ARCH__
        __shared__ float4 cache[DIM];

        for (unsigned int i = 0; i < m_lenses_size; i += DIM) {
            // Load to cache
            const unsigned int idx = i + threadIdx.x;
            if (idx < m_lenses_size) {
                cache[threadIdx.x] = m_lens_int[idx];
            }

            const unsigned int lenses = min(DIM, m_lenses_size - i);
            /*
            if ( blockIdx.x * blockDim.x + threadIdx.x == 0) {
                    printf("\x1B[36mLenses: %d\x1B[0m\n", lenses);
             }
                        */
            __syncthreads();

            // Calculate with cached values
#pragma unroll 32
            for (unsigned int x = 0; x < lenses; x++) {
                float2 t = theta;
                const float4 p1 = cache[x];
                t.x -= p1.z;
                t.y -= p1.w;
                float len = (t.x * t.x) + (t.y * t.y);
                len += p1.x;
                len = __fdividef(p1.y, len);
                t.x *= len;
                t.y *= len;
                alpha.x += t.x;
                alpha.y += t.y;
            }
        }
		if (mass_scale != 0.0) {
			const float scale = m_mass_sheet * mass_scale;
			alpha.x += scale * theta.x;
			alpha.y += scale * theta.y;
		}
#endif
        return alpha;
    }

#else
    inline __host__ __device__ Vector2D<float>
    getAlpha(const Vector2D<float> &theta, const float mass_scale = 0.0) {
        Vector2D<float> alpha(0, 0), movedtheta;
        for (int i = 0; i < m_lenses_size; i++) {
            movedtheta = m_lenses[i].getAlpha(theta);
            alpha += movedtheta;
        }
		if (mass_scale != 0.0) {
			const float scale = m_mass_sheet * mass_scale;
			alpha += (theta * scale);
		}
        return alpha;
    }
#endif
};

/**
 * Builder for CompositeLens.
 */
class CompositeLensBuilder {
  private:
    float m_redshift;
    double m_Dd, m_scale;
    std::vector<Plummer> m_lenses;

  public:
    /**
     * New CompositeLensBuilder.
     */
    CompositeLensBuilder(const float redshift, const double Dd = 0.0,
                         const double scale = 1.0) {
        m_redshift = redshift;
        m_Dd = Dd;
        m_scale = scale;
    }

    float redshift() const { return m_redshift; }

    size_t length() const { return m_lenses.size(); }

    /**
     * Add plummer sublens.
     *
     * @param lens Plummer sublens.
     * @param position Plummer sublens position.
     */
    void addLens(const Plummer &lens);

    /**
     * Get cuda lens. Has internal pointer to cuda memory.
     */
    CompositeLens getCuLens();
    /**
     * Get lens. Has internal pointer to host memory.
     */
    CompositeLens getLens();

    bool operator<(const CompositeLensBuilder &cmp) const {
        return m_redshift < cmp.m_redshift;
    }
};
