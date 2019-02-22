#pragma once

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

    __host__ __device__ Vector2D<T> operator*(const T &scalar) const {
        Vector2D<T> n(m_x * scalar, m_y * scalar);
        return n;
    }
    __host__ __device__ Vector2D<T> &operator*=(const T &scalar) {
		m_x *= scalar;
		m_y *= scalar;
		return *this;
	}
};
