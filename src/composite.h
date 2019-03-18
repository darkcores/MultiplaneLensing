#pragma once

#include "plummer.h"
#include "util/vector2d.h"
#include <vector>

#include <cuda_runtime_api.h>

#include <memory>
#include <utility>

/**
 * Data for lens with position.
 */
class LensData {
  public:
	/**
	 * Plummer lens.
	 */
    Plummer lens;
	/**
	 * Position for the Plummer lens.
	 */
    Vector2D<float> position;
    // bool notlast;
	/**
	 * Create new LensData.
	 *
	 * @param l Lens
	 * @param pos Position.
	 */
    LensData(const Plummer &l, const Vector2D<float> &pos) {
		lens = l;
        position = pos;
        // notlast = true;
    }
	/**
	 * Empty constructor, device only.
	 */
    __host__ __device__ LensData() {
		// notlast = false;
	}
};

/**
 * CompositeLens, contains an array of Plummer lenses.
 */
class CompositeLens {
  private:
    double m_Dd, m_Ds, m_Dds;
    double m_D;
    float m_Df;
    float m_scale;
    LensData *m_data_ptr;
    const LensData *__restrict__ cur_data_ptr;
    int length;

  public:
	/**
	 * New CompositeLens.
	 * 
	 * @param Dd Angular diameter distance.
	 * @param Ds Source angular diameter distance.
	 * @param Dds Lens<->Source angular diameter distance.
	 * @param data_ptr Pointer to LensData.
	 * @param size data_ptr size.
	 * @param scale scale factor single precision.
	 */
    CompositeLens(const double Dd, const double Ds, const double Dds,
                  LensData *data_ptr, size_t size, float scale = 60);
    // __device__ CompositeLens() : cur_data_ptr(nullptr) {}

    __host__ __device__ float distance() const { return m_Dd; }

	/**
	 * Get alpha vector (double precision).
	 *
	 * @param theta Theta vector.
	 * @returns Alpha vector.
	 */
    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;

	/**
	 * Get alpha vector (single precision, with scaling).
	 *
	 * @param theta Theta vector.
	 * @returns Alpha vector.
	 */
    __host__ __device__ Vector2D<float>
    getAlphaf(const Vector2D<float> &theta) const {
        // auto scaledtheta = theta * m_scale;
		// printf("Scale factor %f\n", m_scale);
        Vector2D<float> alpha(0, 0);
        for (int i = 0; i < length; i++) {
            auto movedtheta =
                theta - (cur_data_ptr[i].position);
            alpha += cur_data_ptr[i].lens.getAlphaf(movedtheta * m_scale);
        }
        // Some other tests with less registers but slower (in certain tests)
        // Need to try with the rest to see how this performs
        /*
        auto cur_ptr = cur_data_ptr;
        do {
    auto movedtheta = scaledtheta - (cur_ptr->position * m_scale);
    alpha += cur_ptr->lens.getAlphaf(movedtheta);
                // cur_ptr += sizeof(LensData);
                cur_ptr ++;
        } while(cur_ptr->notlast);
        */
        /*
        int i = 0;
        while (cur_data_ptr[i].notlast) {
    auto movedtheta = scaledtheta - (cur_data_ptr[i].position * m_scale);
    alpha += cur_data_ptr[i].lens.getAlphaf(movedtheta);
                i++;
        }
*/
        // theta /= m_scale;
        alpha /= m_scale;
        return alpha;
    }

	/**
	 * Get beta vector (double precision).
	 *
	 * @param theta Theta vector.
	 * @returns Beta vector.
	 */
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;

	/**
	 * Get beta vector (single precision, with scaling).
	 *  
	 * @param theta Theta vector. 
	 * @returns Beta vector.
	 */
    __host__ __device__ Vector2D<float>
    getBetaf(const Vector2D<float> &theta) const {
        Vector2D<float> beta;
        beta = theta - getAlphaf(theta) * m_Df;
        return beta;
    }

	/**
	 * Get beta vector (single precision, with scaling).
	 *  
	 * @param theta Theta vector. 
	 * @param Ds Source angular distance.
	 * @param Dds lens<->source angular distance.
	 * @returns Beta vector.
	 */
    __host__ __device__ Vector2D<float> getBetaf(const Vector2D<float> &theta,
                                                 const float &Ds,
                                                 const float &Dds) const {
        Vector2D<float> beta;
        beta = theta - getAlphaf(theta) * (Dds / Ds);
        return beta;
    }

	/**
	 * Set masses on sublenses (Host only).
	 *
	 * @param masses Masses for sublenses.
	 */
    __host__ void setMasses(const double *__restrict__ masses) {
        for (int i = 0; i < length; i++) {
            m_data_ptr[i].lens.setMass(masses[i]);
        }
    }

	/**
	 * Set masses on sublenses (Device only).
	 *
	 * @param i Sublens index.
	 * @param mass Mass for sublens.
	 */
	__device__ void setMass(const int i, const double mass) {
		m_data_ptr[i].lens.setMass(mass);
	}

	/*
    __host__ __device__ void setDistance(const double Dd) {
        for (int i = 0; i < length; i++) {
            m_data_ptr[i].lens.setDistance(Dd);
        }
    }

    __host__ __device__ void setSource(const double Ds, const double Dds) {
        m_Ds = Ds;
        m_Dds = Dds;
		m_D = Dds / Ds;
		m_Df = Dds / Ds;
        for (int i = 0; i < length; i++) {
            m_data_ptr[i].lens.setSource(Ds, Dds);
        }
    }
	*/
};

/**
 * Builder for CompositeLens. Please beware currently also manager
 * memory for CompositeLens, will change in the future. For now, do
 * not delete until CompositeLens is no longer needed.
 */
class CompositeLensBuilder {
  private:
	std::vector<LensData> m_lenses;
    LensData *lens_ptr;

    double m_Dd, m_Ds, m_Dds;
	float m_redshift;
    float m_scale;

	bool cuda = false;

  public:
	/**
	 * New CompositeLensBuilder.
	 */
    CompositeLensBuilder() { m_scale = 60; }
	/**
	 * Cleanup cuda allocations if needed.
	 */
	~CompositeLensBuilder() {
		if (cuda)
			cuFree();
	}

	/**
	 * Set lens angular diameter distance.
	 *
	 * @param Dd Angular diameter distance.
	 */
    void setDistance(const double Dd);
	/**
	 * Set redshift for this lens / plane.
	 *
	 * @param z Redshift. 
	 */
	void setRedshift(const float z) { m_redshift = z; }
	/**
	 * Set the source angular distance.
	 *
	 * @param Ds Source angular distance.
	 * @param Dds lens<->source angular distance.
	 */
    void setSource(const double Ds, const double Dds);
	/**
	 * set single precision float scale factor.
	 *
	 * @param scale scale factor.
	 */
    void setScale(const float scale);
	float redshift() const { return m_redshift; }

	/**
	 * Add plummer sublens.
	 *
	 * @param lens Plummer sublens.
	 * @param position Plummer sublens position.
	 */
    void addLens(Plummer &lens, Vector2D<float> position);
    void clear();

	/**
	 * Get cuda lens. Has internal pointer to cuda memory.
	 */
    CompositeLens getCuLens();
	/**
	 * Clear cuda data, automatically called by destructor.
	 */
	void cuFree();
	/**
	 * Get lens. Has internal pointer to host memory.
	 */
    CompositeLens getLens();
};
