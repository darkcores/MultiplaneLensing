#pragma once

#include "util/vector2d.h"
#include <cstdint>
#include <cstring>
#include <vector>

/**
 * Source point data.
 */
class SourceData {
  public:
    /**
     * Source point position.
     */
    Vector2D<float> position;
    /**
     * Source point radius.
     */
    float radius;
    /**
     * Source point color.
     */
    uint8_t color;

    /**
     * New source point data.
     *
     * @param pos Angular position..
     * @param r Radius.
     * @param c Color value (greyscale / map)
     */
    SourceData(const Vector2D<float> pos, const float r,
               const uint8_t c = 255) {
        position = pos;
        radius = r;
        color = c;
    }

    /**
     * Check if a vector is inside this point.
     *
     * @param theta Theta position.
     * @returns color if hit, 0 if none.
     */
    __host__ __device__ uint8_t check_hit(const Vector2D<float> &theta) const {
        float diff = sqrt((theta - position).lengthSq());
        diff = fabs(diff);
#ifdef __CUDACC__
// printf("Point diff: %f\n", diff);
#endif
        if (diff < radius) {
            return color;
        }
        return 0;
    }
};

/**
 * Source plane class.
 */
class SourcePlane {
  private:
    float m_redshift, m_Ds, m_Dds;
    const SourceData *__restrict__ m_points;
    SourceData *m_points_ptr;
    int m_points_length;
    bool m_cuda;

  public:
    /**
     * New SourcePlane.
     *
     * @param redshift Redshift distance.
     * @param points SourceData points.
     * @param points_length Size of SourceData points.
     */
    SourcePlane(const float redshift, SourceData *points,
                const int points_length, bool cuda = false)
        : m_points(points) {
        m_points_length = points_length;
        m_redshift = redshift;
        m_points_ptr = points;
        m_cuda = cuda;
    }

    /**
     * Memory cleanup.
     */
    int destroy() {
        if (m_points_ptr) {
            if (m_cuda) {
#ifdef __CUDACC__
                cudaFree(m_points_ptr);
#endif
            } else {
                free(m_points_ptr);
            }
            m_points_ptr = NULL;
        }
        return 0;
    }

    /**
     * Set the source angular distance.
     *
     * @param Ds Source angular distance.
     * @param Dds lens<->source angular distance.
     */
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

    /**
     * Check if theta hits a source point in this plane.
     *
     * @param theta Theta vector.
     */
    __host__ __device__ uint8_t check_hit(const Vector2D<float> &theta) const {
        // Check for each point |diff| < radius
#ifdef __CUDACC__
// printf("Points in plane: %d\n", m_points_length);
#endif
        for (int i = 0; i < m_points_length; i++) {
            uint8_t p = m_points[i].check_hit(theta);
            if (p > 0)
                return p;
        }
        return 0;
    }
};

/**
 * Builder class for SourcePlane.
 */
class SourcePlaneBuilder {
  private:
    float m_redshift;
    std::vector<SourceData> m_points;
    bool cuda = false;
    SourceData *ptr;

  public:
    /**
     * New SourcePlaneBuilder.
     *
     * @param redshift Redshift distance.
     */
    SourcePlaneBuilder(const float redshift) { m_redshift = redshift; }

    /**
     * Add source point data.
     *
     * @param point Angular position..
     * @param radius Radius.
     * @param color Color value (greyscale / map).
     */
    void addPoint(const Vector2D<float> point, const float radius,
                  const uint8_t color = 255) {
        m_points.push_back(SourceData(point, radius, color));
    }

    /**
     * Get LensPlane for CPU.
     */
    SourcePlane getPlane() {
        ptr = (SourceData *)malloc(sizeof(SourceData) * m_points.size());
		memcpy(ptr, &m_points[0], sizeof(SourceData) * m_points.size());

        return SourcePlane(m_redshift, ptr, m_points.size());
    }

    /**
     * Get LensPlane for GPU/CUDA.
     */
    SourcePlane getCuPlane() {
        cuda = true;
#ifdef __CUDACC__
        // dev_m_points = m_points;
        // ptr = thrust::raw_pointer_cast(&dev_m_points[0]);
        cudaMalloc(&ptr, sizeof(SourceData) * m_points.size());
        cudaMemcpy(ptr, &m_points[0], sizeof(SourceData) * m_points.size(),
                   cudaMemcpyHostToDevice);
#else
        ptr = nullptr;
#endif
        return SourcePlane(m_redshift, ptr, m_points.size(), true);
    }
};
