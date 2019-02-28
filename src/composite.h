#pragma once

#include "plummer.h"
#include "util/vector2d.h"

#include <cuda_runtime_api.h>
#ifdef __CUDACC__
#include <thrust/device_vector.h>
#endif
#include <thrust/host_vector.h>

#include <memory>
#include <utility>

class LensData {
  public:
    Plummer lens;
    Vector2D<float> position;
    LensData(const Plummer &l, const Vector2D<float> &pos) : lens(l) {
        position = pos;
    }
};

class CompositeLens {
  private:
    double m_Dd, m_Ds, m_Dds;
	double m_D;
	float m_Df;
    float m_scale;
    LensData *cur_data_ptr;
    size_t length;

  public:
    CompositeLens(const double Dd, const double Ds, const double Dds,
                  LensData *data_ptr, size_t size, float scale = 60);

    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float> getAlphaf(Vector2D<float> theta) const;
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float> getBetaf(Vector2D<float> theta) const;
};

class CompositeLensBuilder {
  private:
#ifdef __CUDACC__
    thrust::device_vector<LensData> dev_m_lenses;
#endif
    thrust::host_vector<LensData> m_lenses;

    double m_Dd, m_Ds, m_Dds;
    float m_scale;

  public:
    CompositeLensBuilder();

    void setDistance(const double Dd);
    void setSource(const double Ds, const double Dds);
    void setScale(const float scale = 60);

    void addLens(Plummer &lens, Vector2D<float> position);
    void prepare();
    void clear();

    CompositeLens getCuLens();
    CompositeLens getLens();
};
