#include "gtest/gtest.h"

#include <composite.h>
#include <multiplane.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <fstream>
#include <thrust/device_vector.h>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass) {
    CompositeLensBuilder lensbuilder;
    double xstart = -width / 2;
    double xend = width / 2;
    double xstep = width / (N - 1);
    double ystart = -height / 2;
    double yend = height / 2;
    double ystep = height / (N - 1);
    // This works for this test, but might not always work TODO
    for (double x = xstart; x <= xend; x += xstep) {
        for (double y = ystart; y <= yend; y += ystep) {
            Plummer plum(Dd, mass, angularwidth);
            lensbuilder.addLens(plum, Vector2D<float>(x, y));
        }
    }
    return lensbuilder;
}

std::vector<Vector2D<float>> thetaGrid(Vector2D<float> topleft,
                                       Vector2D<float> botright,
                                       int numX = 1000, int numY = 1000) {
    std::vector<Vector2D<float>> points;
    Vector2D<float> step = botright - topleft;
    step.setX(fabs(step.x() / numX));
    step.setY(fabs(step.y() / numY));
    float y = topleft.y();
    for (int iy = 0; iy < numY; iy++) {
        float x = topleft.x();
        for (int ix = 0; ix < numX; ix++) {
            points.push_back(Vector2D<float>(x, y));
            x += step.x();
        }
        y -= step.y();
    }
    return points;
}

__global__ void traceThetaKernel(const int n, const Multiplane mp,
                                 const float *__restrict__ xpoints,
                                 const float *__restrict__ ypoints,
                                 uint8_t *__restrict__ output) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        Vector2D<float> vec(xpoints[i], ypoints[i]);
        uint8_t p = mp.traceTheta(vec);
        output[i] = p;
    }
}

TEST(CuMultiplaneTests, TestBetafUberKernel) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
    double z_d2 = 0.8;
    double z_s = 1.2;
    double z_s2 = 0.6;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    auto lensbuilder = createGrid(Dd1, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
    lensbuilder.setRedshift(z_d1);
	lensbuilder.setScale(10);
    auto lensbuilder2 = createGrid(Dd2, 6, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                   5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
    lensbuilder2.setRedshift(z_d2);
	lensbuilder2.setScale(10);

    MultiplaneBuilder planebuilder(cosm);
    planebuilder.addPlane(&lensbuilder);
    planebuilder.addPlane(&lensbuilder2);

    SourcePlaneBuilder sourcebuilder(z_s);
    sourcebuilder.addPoint(Vector2D<float>(1 * ANGLE_ARCSEC, -2 * ANGLE_ARCSEC),
                           1 * ANGLE_ARCSEC);
    sourcebuilder.addPoint(Vector2D<float>(-3 * ANGLE_ARCSEC, 2 * ANGLE_ARCSEC),
                           1 * ANGLE_ARCSEC);
    sourcebuilder.addPoint(Vector2D<float>(1 * ANGLE_ARCSEC, -2 * ANGLE_ARCSEC),
                           1 * ANGLE_ARCSEC);
    auto sourceplane = sourcebuilder.getCuPlane();
    planebuilder.addSourcePlane(sourceplane);

    SourcePlaneBuilder sourcebuilder2(z_s2);
    sourcebuilder2.addPoint(
        Vector2D<float>(-4 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC), 1 * ANGLE_ARCSEC,
        128);
    sourcebuilder2.addPoint(
        Vector2D<float>(-3 * ANGLE_ARCSEC, -10 * ANGLE_ARCSEC),
        1 * ANGLE_ARCSEC, 128);
    sourcebuilder2.addPoint(
        Vector2D<float>(-3 * ANGLE_ARCSEC, -4 * ANGLE_ARCSEC), 1 * ANGLE_ARCSEC,
        128);
    auto sourceplane2 = sourcebuilder2.getCuPlane();
    planebuilder.addSourcePlane(sourceplane2);

    auto multiplane = planebuilder.getCuMultiPlane();
    // Multiplane *mp;
    // cudaMalloc(&mp, sizeof(Multiplane));
    // cudaMemcpy(mp, &multiplane, sizeof(Multiplane), cudaMemcpyHostToDevice);

    // Points and layout for cuda
    auto points =
        thetaGrid(Vector2D<float>(-30 * ANGLE_ARCSEC, 30 * ANGLE_ARCSEC),
                  Vector2D<float>(30 * ANGLE_ARCSEC, -30 * ANGLE_ARCSEC));
    thrust::host_vector<float> xpoints;
    thrust::host_vector<float> ypoints;
    for (auto &p : points) {
        xpoints.push_back(p.x());
        ypoints.push_back(p.y());
    }
    thrust::device_vector<float> dev_x(xpoints);
    thrust::device_vector<float> dev_y(ypoints);
    float *dev_x_ptr = thrust::raw_pointer_cast(&dev_x[0]);
    float *dev_y_ptr = thrust::raw_pointer_cast(&dev_y[0]);

    auto size = xpoints.size();
    thrust::device_vector<uint8_t> dev_o(size);
    uint8_t *dev_o_ptr = thrust::raw_pointer_cast(&dev_o[0]);

    std::cout << "Setup done, starting kernel" << std::endl;
    // size_t newHeapSize = 1024 * 1000 * 32;
    // cudaDeviceSetLimit(cudaLimitMallocHeapSize, newHeapSize);
    traceThetaKernel<<<(size / 64) + 1, 64>>>(size, multiplane, dev_x_ptr, dev_y_ptr,
                                        dev_o_ptr);

    std::cout << "Kernel: " << cudaGetErrorString(cudaPeekAtLastError())
              << std::endl;
    cudaDeviceSynchronize();
    std::cout << "Kernel done" << std::endl;
    thrust::host_vector<uint8_t> output = dev_o;

    std::ofstream testimg("testimage_cu.raw", std::ios::binary);
    testimg.write((char *)&output[0], sizeof(uint8_t) * output.size());
    testimg.close();
}

TEST(CuMultiplaneTests, TestBetafMultiKernel) {}