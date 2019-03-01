#include "gtest/gtest.h"

#include <composite.h>
#include <thrust/device_vector.h>
#include <util/constants.h>
#include <util/cosmology.h>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass, double Ds = 0,
                                double Dds = 0) {
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
    lensbuilder.setSource(Ds, Dds);
    return lensbuilder;
}

__global__ void alphatest(int n, CompositeLens *lens,
                          Vector2D<double> *alphas) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    Vector2D<double> vec((i % 7) * ANGLE_ARCSEC, (i % 5) * ANGLE_ARCSEC);
    alphas[i] = lens->getAlpha(vec);
}

__global__ void betatest(int n, CompositeLens *lens, Vector2D<double> *betas) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    Vector2D<double> vec((i % 7) * ANGLE_ARCSEC, (i % 5) * ANGLE_ARCSEC);
    betas[i] = lens->getBeta(vec);
}

__global__ void betaftest(int n, const CompositeLens *const lens, float *thetax,
                          float *thetay, float *betax, float *betay) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ CompositeLens slens;
	if (threadIdx.x == 0) {
		slens = *lens;
	}
	__syncthreads();
    // Vector2D<float> vec((i % 7) * ANGLE_ARCSEC, (i % 5) * ANGLE_ARCSEC);
	Vector2D<float> vec(thetax[i], thetay[i]);
    auto betas = slens.getBetaf(vec);
	betax[i] = betas.x();
	betay[i] = betas.y();
}

TEST(CompositeCuTests, TestAlpha) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto lensbuilder = createGrid(Dd, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
    // lens.prepare();
    thrust::device_vector<CompositeLens> dv;
    dv.push_back(lensbuilder.getCuLens());
    auto l_ptr = thrust::raw_pointer_cast(&dv[0]);

    thrust::device_vector<Vector2D<double>> alphas(1048576);
    auto a_ptr = thrust::raw_pointer_cast(&alphas[0]);

    alphatest<<<1048576 / 256, 256>>>(1048576, l_ptr, a_ptr);
    thrust::host_vector<Vector2D<double>> h_alphas(alphas);

    /*
    auto tlens = lensbuilder.getLens();
for (int i = 0; i < 128; i += 64) {
    for (int j = 0; j < 64; j++) {
        int idx = i + j;
        Vector2D<double> vec((idx % 7) * ANGLE_ARCSEC,
                             (idx % 5) * ANGLE_ARCSEC);
        Vector2D<double> result = tlens.getAlpha(vec);
        EXPECT_EQ(h_alphas[idx], result);
    }
}
    */
}

TEST(CompositeCuTests, TestBeta) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    double z_s = 1.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    auto lensbuilder = createGrid(Dd, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR, Ds, Dds);
    auto lens = lensbuilder.getCuLens();

    thrust::device_vector<CompositeLens> dv;
    dv.push_back(lensbuilder.getCuLens());
    auto l_ptr = thrust::raw_pointer_cast(&dv[0]);

    thrust::device_vector<Vector2D<double>> betas(1048576);
    auto b_ptr = thrust::raw_pointer_cast(&betas[0]);

    betatest<<<1048576 / 256, 256>>>(1048576, l_ptr, b_ptr);
    thrust::host_vector<Vector2D<double>> h_betas(betas);
}

TEST(CompositeCuTests, TestBetaF) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    double z_s = 1.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    auto lensbuilder = createGrid(Dd, 10, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR, Ds, Dds);
    auto lens = lensbuilder.getCuLens();

    thrust::device_vector<CompositeLens> dv;
    dv.push_back(lensbuilder.getCuLens());
    auto l_ptr = thrust::raw_pointer_cast(&dv[0]);

    thrust::host_vector<float> thetax;
    thrust::host_vector<float> thetay;
    for (size_t i = 0; i < 16777216; i++) {
        thetax.push_back((i % 7) * ANGLE_ARCSEC);
        thetay.push_back((i % 5) * ANGLE_ARCSEC);
    }
    thrust::device_vector<float> dev_thetax(thetax);
    thrust::device_vector<float> dev_thetay(thetay);
    thrust::device_vector<float> betax(16777216);
    thrust::device_vector<float> betay(16777216);
    auto bx_ptr = thrust::raw_pointer_cast(&betax[0]);
    auto by_ptr = thrust::raw_pointer_cast(&betay[0]);
    auto tx_ptr = thrust::raw_pointer_cast(&dev_thetax[0]);
    auto ty_ptr = thrust::raw_pointer_cast(&dev_thetay[0]);

	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    betaftest<<<16777216 / 512, 512>>>(16777216, l_ptr, tx_ptr, tx_ptr, bx_ptr, by_ptr);
    // thrust::host_vector<Vector2D<float>> h_betas(betas);
}
