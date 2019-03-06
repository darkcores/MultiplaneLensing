#pragma once

#include "util/vector2d.h"
#include <cstdint>
#include <vector>

#include <cuda_runtime_api.h>

class SourceData {
  public:
    Vector2D<float> position;
    float radius;
    uint8_t color;

    SourceData(const Vector2D<float> pos, const float r,
               const uint8_t c = 255) {
        position = pos;
        radius = r;
        color = c;
    }
};

class SourcePlane {
  private:
    float m_redshift, m_Ds, m_Dds;
    const SourceData *__restrict__ m_points;
    int m_points_length;

  public:
    SourcePlane(const float redshift, const SourceData *points,
                const int points_length)
        : m_points(points) {
        m_points_length = points_length;
        m_redshift = redshift;
    }

    void setSource(float Ds, float Dds) {
        m_Ds = Ds;
        m_Dds = Dds;
    }

    __host__ __device__ const float &redshift() const { return m_redshift; }
    __host__ __device__ const float &ds() const { return m_Ds; }
    __host__ __device__ const float &dds() const { return m_Dds; }

    __host__ __device__ bool operator<(const SourcePlane &cmp) const {
        return m_redshift < cmp.redshift();
    }

    __host__ __device__ uint8_t check_hit(const Vector2D<float> &theta) const;
};

class SourcePlaneBuilder {
  private:
    float m_redshift;
    std::vector<SourceData> m_points;
    SourceData *ptr;

  public:
    SourcePlaneBuilder(const float redshift) { m_redshift = redshift; }

    void addPoint(const Vector2D<float> point, const float radius) {
        m_points.push_back(SourceData(point, radius));
    }

    SourcePlane getPlane() const {
        return SourcePlane(m_redshift, &m_points[0], m_points.size());
    }

    SourcePlane getCuPlane() {
#ifdef __CUDACC__
        // dev_m_points = m_points;
        // ptr = thrust::raw_pointer_cast(&dev_m_points[0]);
        cudaMalloc(&ptr, sizeof(SourceData) * m_points.size());
        cudaMemcpy(ptr, &m_points[0], sizeof(SourceData) * m_points.size(),
                   cudaMemcpyHostToDevice);
#else
        ptr = nullptr;
#endif
        return SourcePlane(m_redshift, ptr, m_points.size());
    }

    void cuFree() {
#ifdef __CUDACC__
        cudaFree(ptr);
#endif
    }
};
