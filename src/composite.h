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

  public:
    CompositeLens(Plummer *lenses, const int lenses_size,
                  const bool cuda = false);

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
    __host__ __device__ float2 getAlpha(const float2 &theta) const {
		const unsigned int DIM = 128;
        float2 alpha, t;
        alpha.x = 0;
        alpha.y = 0;
#ifdef __CUDA_ARCH__
        __shared__ float4 cache[DIM];

        // #pragma unroll 64
        for (unsigned int i = 0; i < m_lenses_size; i += DIM) {
            // Load to cache
            const unsigned int idx = i + threadIdx.x;
            if (idx < m_lenses_size) {
                // const float4 *r = m_lens_int + (i + threadIdx.x);
                cache[threadIdx.x] = m_lens_int[idx];
                /*
asm("ld.global.ca.v4.f32 {%0, %1, %2, %3}, [%4];"
: "=f"(p1.x), "=f"(p1.y), "=f"(p1.z), "=f"(p1.w)
: "l"(r));
                */
            }

            const unsigned int lenses = min(DIM, m_lenses_size - i);
            // Calculate with cached values
            // p1 = m_lens_int[i];
            /*
            if ( blockIdx.x * blockDim.x + threadIdx.x == 0) {
                    printf("\x1B[36mLenses: %d\x1B[0m\n", lenses);
            }
			*/
            __syncthreads();
#pragma unroll 32
            for (unsigned int x = 0; x < lenses; x++) {
                t = theta;
                const float4 p1 = cache[x];
                t.x -= p1.z;
                t.y -= p1.w;
                float len = (t.x * t.x) + (t.y * t.y);
                len += p1.x;
                // len = (1 / len);
				len = __fdividef(1.0, len);
                len *= p1.y;
                t.x *= len;
                t.y *= len;
                alpha.x += t.x;
                alpha.y += t.y;
            }
        }
#endif
        return alpha;
    }

#else
    inline __host__ __device__ Vector2D<float> getAlpha(const Vector2D<float> &theta) {
        Vector2D<float> alpha(0, 0), movedtheta;
        for (int i = 0; i < m_lenses_size; i++) {
            movedtheta = m_lenses[i].getAlpha(theta);
            alpha += movedtheta;
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
    std::vector<Plummer> m_lenses;

  public:
    /**
     * New CompositeLensBuilder.
     */
    CompositeLensBuilder(const float redshift) { m_redshift = redshift; }

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
