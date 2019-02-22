#include "gtest/gtest.h"

#include <util/vector2d.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>

__global__ void scalegpu(int n, Vector2D<double> *data, double scalar) {
    int i = threadIdx.x;

    if (i < n) {
        data[i].setX(i);
        data[i].setY(i);
        data[i] *= scalar;
    }
}

template <typename T> struct mul {
    __host__ __device__ T operator()(const T &x, double y) const {
        return x * y;
    }
};

TEST(Vector2DCuTests, ScalingDouble) {
    thrust::host_vector<Vector2D<double>> H(4);
    thrust::device_vector<Vector2D<double>> D = H;
    auto raw = thrust::raw_pointer_cast(&D[0]);
    scalegpu<<<1, 32>>>(4, raw, 4.0);
    mul<Vector2D<double>> op;
    thrust::transform(D.begin(), D.end(), thrust::make_constant_iterator(0.25),
                      D.begin(), op);
    thrust::copy(D.begin(), D.end(), H.begin());
	ASSERT_EQ(H[0].x(), 0);
	ASSERT_EQ(H[1].x(), 1);
}