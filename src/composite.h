#pragma once

#include "plummer.h"
#include <vector>

class CompositeLens {
  private:
    Plummer *__restrict__ m_lenses;
    const int m_lenses_size;
    const bool m_cuda;

  public:
    CompositeLens(Plummer *lenses, const int lenses_size,
                  const bool cuda = false)
        : m_lenses(lenses), m_lenses_size(lenses_size), m_cuda(cuda) {}

    int destroy();

	__host__ void update(const float *__restrict__ factors) {
		for (int i = 0; i < m_lenses_size; i++) {
			m_lenses[i].update(factors[i]);
		}
	}

	__device__ void update(const float &factor, const int idx) {
		m_lenses[idx].update(factor);
	}
	
#ifdef __CUDACC__
    /**
     * Get alpha vector (single precision, with scaling).
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlpha(const float2 &theta) const {
        float2 alpha, movedtheta;
        alpha.x = 0;
        alpha.y = 0;
#ifdef __CUDA_ARCH__
#pragma unroll 16
#endif
        for (int i = 0; i < m_lenses_size; i++) {
            movedtheta = m_lenses[i].getAlpha(theta);
            alpha.x += movedtheta.x;
            alpha.y += movedtheta.y;
        }
        return alpha;
    }
#else
	__host__ __device__ Vector2D<float> getAlpha(const Vector2D<float> &theta) {
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
 * Builder for CompositeLens. Please beware currently also manager
 * memory for CompositeLens, will change in the future. For now, do
 * not delete until CompositeLens is no longer needed.
 */
class CompositeLensBuilder {
  private:
    std::vector<Plummer> m_lenses;

    float m_redshift;

  public:
    /**
     * New CompositeLensBuilder.
     */
    CompositeLensBuilder(const float redshift) {
        m_redshift = redshift;
    }

    float redshift() const { return m_redshift; }

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
