#pragma once

#include <cmath>
// CUDA includes
#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#else
#ifndef __host__
#define __host__
#define __device__
#endif
#endif

/**
 * Vector2D template class.
 */
template <typename T> class Vector2D {
  private:
    T m_x, m_y;

  public:
    /**
     * Create a new empty Vector2D object.
     */
    __host__ __device__ Vector2D() {}
    /**
     * Create a new empty Vector2D object.
     *
     * @param x X component.
     * @param y Y component.
     */
    __host__ __device__ Vector2D(const T &__restrict__ x,
                                 const T &__restrict__ y) {
        m_x = x;
        m_y = y;
    };

    __host__ __device__ inline const T &x() const { return m_x; }
    __host__ __device__ inline const T &y() const { return m_y; }

    __host__ __device__ inline void setX(const T &x) { m_x = x; }
    __host__ __device__ inline void setY(const T &y) { m_y = y; }

    __host__ __device__ inline Vector2D<T> operator*(const T &scalar) const {
        return Vector2D<T>(m_x * scalar, m_y * scalar);
    }
    __host__ __device__ inline void operator*=(const T &scalar) {
        m_x *= scalar;
        m_y *= scalar;
        // return *this;
    }
    __host__ __device__ inline void operator+=(const Vector2D<T> &add) {
        m_x += add.x();
        m_y += add.y();
        // return *this;
    }

    __host__ __device__ inline Vector2D<T> operator/(const T &scalar) const {
        return Vector2D<T>(m_x / scalar, m_y / scalar);
    }
    __host__ __device__ inline void operator/=(const T &scalar) {
        m_x /= scalar;
        m_y /= scalar;
        // return *this;
    }

    template <typename Y>
    __host__ __device__ inline Vector2D<T>
    operator-(const Vector2D<Y> &diff) const {
        return Vector2D<T>(m_x - diff.x(), m_y - diff.y());
    }
    __host__ __device__ inline void operator-=(const Vector2D<T> &diff) {
        m_x -= diff.x();
        m_y -= diff.y();
        // return *this;
    }
#ifdef __CUDACC__
    __host__ __device__ inline void operator-=(const float2 &diff) {
        m_x -= diff.x;
        m_y -= diff.y;
        // return *this;
    }
#endif
    __host__ __device__ inline void rdiff(const Vector2D<T> &diff) {
        m_x = diff.x() - m_x;
        m_y = diff.y() - m_y;
    }

    __host__ __device__ inline bool operator==(const Vector2D<T> &cmp) const {
        return (m_x == cmp.x() && m_y == cmp.y());
    }
    __host__ __device__ Vector2D<T> pow(const int exp) const {
        Vector2D<T> x(pow(m_x, exp), pow(m_y, exp));
        return x;
    }

    /**
     * Get Length squared.
     *
     * @returns length.
     */
    __host__ __device__ inline T lengthSq() const {
        return m_x * m_x + m_y * m_y;
    }
};
