#pragma once

#include <cmath>
// CUDA includes
#include <cuda_runtime_api.h>

template <typename T> class Vector2D {
  private:
    T m_x, m_y;

  public:
    __host__ __device__ Vector2D() {}
    __host__ __device__ Vector2D(T x, T y) {
        m_x = x;
        m_y = y;
    };

    __host__ __device__ const T &x() const { return m_x; }
    __host__ __device__ const T &y() const { return m_y; }

    __host__ __device__ void setX(const T &x) { m_x = x; }
    __host__ __device__ void setY(const T &y) { m_y = y; }

    template <typename Y>
    __host__ __device__ Vector2D<T> operator*(const Y &scalar) const {
        Vector2D<T> n(m_x * scalar, m_y * scalar);
        return n;
    }
    template <typename Y>
    __host__ __device__ Vector2D<T> &operator*=(const Y &scalar) {
        m_x *= scalar;
        m_y *= scalar;
        return *this;
    }

    template <typename Y>
    __host__ __device__ Vector2D<T> operator/(const Y &scalar) const {
        Vector2D<T> n(m_x / scalar, m_y / scalar);
        return n;
    }

    __host__ __device__ Vector2D<T> operator-(const Vector2D<T> &diff) const {
        Vector2D<T> n(m_x - diff.x(), m_y - diff.y());
        return n;
    }
    __host__ __device__ Vector2D<T> &operator-=(const Vector2D<T> &diff) {
        m_x -= diff.x();
        m_y -= diff.y();
        return *this;
    }

    __host__ __device__ T length() const {
        T absX = abs(m_x);
        T absY = abs(m_y);

        if (absX > absY) {
            T tmp = absY / absX;
            return absX * std::sqrt((T)1.0 + tmp * tmp);
        }
        // absX <= absY
        if (absY == 0) // => absx == 0
            return 0;
        T tmp = absX / absY;
        return absY * std::sqrt(tmp * tmp + (T)1.0);
        // return sqrt(m_x * m_x + m_y * m_y);
    }

    __host__ __device__ T lengthSq() const { return m_x * m_x + m_y * m_y; }
};
