#pragma once

#include "plummer.h"
#include "util/vector2d.h"

#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
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
    thrust::host_vector<LensData> m_lenses;
    thrust::device_vector<LensData> dev_m_lenses;
	LensData *cur_data_ptr;
	size_t length;

  public:
    CompositeLens();

    void addLens(Plummer &lens, Vector2D<float> position);
	void prepare();
    void clear();

    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float> getAlphaf(Vector2D<float> theta) const;
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float> getBetaf(Vector2D<float> theta) const;

    __host__ __device__ void setDistance(const double Dd);
    __host__ __device__ void setSource(const double Ds, const double Dds);
    __host__ __device__ void setScale(const float scale = 60);
};
