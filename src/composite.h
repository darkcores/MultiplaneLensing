#pragma once

#include "plummer.h"
#include "util/vector2d.h"
#include <vector>

#include <cuda_runtime_api.h>

#include <memory>
#include <utility>

class LensData {
  public:
    Plummer lens;
    Vector2D<float> position;
    // bool notlast;
    LensData(const Plummer &l, const Vector2D<float> &pos) {
		lens = l;
        position = pos;
        // notlast = true;
    }
    __host__ __device__ LensData() {
		// notlast = false;
	}
};

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
    CompositeLens(const double Dd, const double Ds, const double Dds,
                  LensData *data_ptr, size_t size, float scale = 60);
    // __device__ CompositeLens() : cur_data_ptr(nullptr) {}

    __host__ __device__ float distance() const { return m_Dd; }

    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;
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
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float>
    getBetaf(const Vector2D<float> &theta) const {
        Vector2D<float> beta;
        beta = theta - getAlphaf(theta) * m_Df;
        return beta;
    }

    __host__ __device__ Vector2D<float> getBetaf(const Vector2D<float> &theta,
                                                 const float &Ds,
                                                 const float &Dds) const {
        Vector2D<float> beta;
        beta = theta - getAlphaf(theta) * (Dds / Ds);
        return beta;
    }

    __host__ __device__ void setMasses(const double *__restrict__ masses) {
        for (int i = 0; i < length; i++) {
            m_data_ptr[i].lens.setMass(masses[i]);
        }
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

class CompositeLensBuilder {
  private:
	std::vector<LensData> m_lenses;
    LensData *lens_ptr;

    double m_Dd, m_Ds, m_Dds;
	float m_redshift;
    float m_scale;

	bool cuda = false;

  public:
    CompositeLensBuilder() { m_scale = 60; }
	~CompositeLensBuilder() {
		if (cuda)
			cuFree();
	}

    void setDistance(const double Dd);
	void setRedshift(const float z) { m_redshift = z; }
    void setSource(const double Ds, const double Dds);
    void setScale(const float scale);
	float redshift() const { return m_redshift; }

    void addLens(Plummer &lens, Vector2D<float> position);
    void clear();

    CompositeLens getCuLens();
	void cuFree();
    CompositeLens getLens();
};
