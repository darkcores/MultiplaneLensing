#pragma once

#include "composite.h"

#include <cuda_runtime_api.h>
#ifdef __CUDACC__
#include <thrust/device_vector.h>
#endif
#include <thrust/host_vector.h>

class PlaneData {
  public:
    PlaneData(CompositeLens lens, double z);
    __device__ PlaneData() {}

    CompositeLens lens;
    float redshift;

    bool operator<(const PlaneData &cmp) const {
        return redshift < cmp.redshift;
    }
};

class SourceData {
  public:
    Vector2D<float> position;
    float radius;

    SourceData(const Vector2D<float> pos, const float r) {
        position = pos;
        radius = r;
    }
};

class SourcePlane {
  private:
    float m_redshift;
    const SourceData *__restrict__ m_points;

  public:
    SourcePlane(const float redshift, const SourceData *points)
        : m_points(points) {
        m_redshift = redshift;
    }

    __host__ __device__ const float &redshift() const { return m_redshift; }

    __host__ __device__ bool operator<(const SourcePlane &cmp) const {
        return m_redshift < cmp.redshift();
    }
};

class SourcePlaneBuilder {
  private:
    float m_redshift;
    thrust::host_vector<SourceData> m_points;
#ifdef __CUDACC__
    thrust::device_vector<SourceData> dev_m_points;
#endif

  public:
    SourcePlaneBuilder(const float redshift) { m_redshift = redshift; }

    void addPoint(const Vector2D<float> point, const float radius) {
        m_points.push_back(SourceData(point, radius));
    }
};

class Multiplane {
  private:
    const PlaneData *__restrict__ m_plane_ptr;
    const SourceData *__restrict__ m_src_plane_ptr;
    int m_plane_length;
    int m_src_length;

  public:
    Multiplane();

  private:
    __host__ __device__ Vector2D<double> getBeta(const Vector2D<double> &theta,
                                                 const int &lens) const;
};

class MultiplaneBuilder {
  private:
    thrust::host_vector<PlaneData> m_data;
    thrust::host_vector<SourcePlane> m_src_data;
#ifdef __CUDACC__
    thrust::device_vector<PlaneData> dev_m_data;
    thrust::device_vector<SourcePlane> dev_m_src_data;
#endif

    void prepare();

  public:
    MultiplaneBuilder();

    void addPlane(const CompositeLens &lens, double z);
    void addSourcePlane(const SourcePlane &plane);

    Multiplane getMultiPlane();
    Multiplane getCuMultiPlane();
};
