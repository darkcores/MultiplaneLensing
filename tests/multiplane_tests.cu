#include "gtest/gtest.h"

#include <fstream>
#include <multiplane.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

CompositeLensBuilder createGrid(float z, double Dd, int N, double width,
                                double height, double angularwidth,
                                double mass) {
    CompositeLensBuilder lensbuilder(z);
    double xstart = -width / 2;
    double xend = width / 2;
    double xstep = width / (N - 1);
    double ystart = -height / 2;
    double yend = height / 2;
    double ystep = height / (N - 1);
    // This works for this test, but might not always work TODO
    for (double x = xstart; x <= xend; x += xstep) {
        for (double y = ystart; y <= yend; y += ystep) {
            Plummer plum(Dd, mass, angularwidth, 1 / ANGLE_ARCSEC,
                         float2{.x = (float)x, .y = (float)y});
            lensbuilder.addLens(plum);
        }
    }
    return lensbuilder;
}

TEST(MultiplaneTests, TestTrace) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
    double z_d2 = 1.0;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    MultiplaneBuilder builder(cosm);
    auto lensbuilder = createGrid(z_d1, Dd1, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder);
    auto lensbuilder2 = createGrid(z_d2, Dd2, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder2);

    std::vector<float> srcs{0.2, 0.8, 1.5, 2.5};
    builder.setRedshifts(srcs);

    auto mp = builder.getCuMultiPlane();

    float2 point{.x = 2, .y = 1};
    thrust::host_vector<float2> thetas;
    for (int i = 0; i < 340; i++) {
        thetas.push_back(point);
    }
    thrust::host_vector<float2> betas(340);

    thrust::device_vector<float2> d_thetas(thetas);
    thrust::device_vector<float2> d_betas(340);

    float2 *thetasptr = thrust::raw_pointer_cast(&d_thetas[0]);
    float2 *betasptr = thrust::raw_pointer_cast(&d_betas[0]);

    mp.traceThetas(thetasptr, betasptr, 340, 0);
    betas = d_betas;
    for (int i = 0; i < 340; i++) {
        ASSERT_EQ(betas[i].x, 2.0);
        ASSERT_EQ(betas[i].y, 1.0);
    }
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 340, 1);
    betas = d_betas;
    for (int i = 0; i < 340; i++) {
        ASSERT_LT(fabs(betas[i].x + 1.579611), 1e-5);
        ASSERT_LT(fabs(betas[i].y + 0.840784), 1e-5);
    }
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 340, 2);
    betas = d_betas;
    for (int i = 0; i < 340; i++) {
        ASSERT_LT(fabs(betas[i].x + 1.736920), 1e-5);
        ASSERT_LT(fabs(betas[i].y + 0.876505), 1e-5);
    }
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 340, 3);
    betas = d_betas;
    for (int i = 0; i < 340; i++) {
        ASSERT_LT(fabs(betas[i].x + 1.310565), 1e-5);
        ASSERT_LT(fabs(betas[i].y + 0.621893), 1e-5);
    }
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.destroy();
}

TEST(MultiplaneTests, TestTraceUpdate) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
    double z_d2 = 1.0;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    MultiplaneBuilder builder(cosm);
    auto lensbuilder = createGrid(z_d1, Dd1, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder);
    auto lensbuilder2 = createGrid(z_d2, Dd2, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder2);

    std::vector<float> srcs{0.2, 0.8, 1.5, 2.5};
    builder.setRedshifts(srcs);

    auto mp = builder.getCuMultiPlane();

    std::vector<std::vector<float>> factors{{0, 1, 0, 1, 0, 0.5, 0, 2, 0},
                                            {1, 0.5, 0, 1, 1.5, 0, 0, 1, 0.5}};
    mp.updateMassesCu(factors);

    float2 point{.x = 2, .y = 1};
    thrust::host_vector<float2> thetas;
    thetas.push_back(point);
    thrust::host_vector<float2> betas(1);

    thrust::device_vector<float2> d_thetas(thetas);
    thrust::device_vector<float2> d_betas(1);

    float2 *thetasptr = thrust::raw_pointer_cast(&d_thetas[0]);
    float2 *betasptr = thrust::raw_pointer_cast(&d_betas[0]);

    mp.traceThetas(thetasptr, betasptr, 1, 0);
    betas = d_betas;
    EXPECT_EQ(betas[0].x, 2.0);
    EXPECT_EQ(betas[0].y, 1.0);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 1);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x - 4.607982), 1e-5);
    EXPECT_LT(fabs(betas[0].y + 1.687718), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 2);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x - 3.180912), 1e-5);
    EXPECT_LT(fabs(betas[0].y + 2.014467), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 3);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x - 1.663422), 1e-5);
    EXPECT_LT(fabs(betas[0].y + 1.857666), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.destroy();
}

TEST(MultiplaneTests, TestTraceFileLens) {
    const double unit = ANGLE_ARCSEC;
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
    double z_d2 = 1.0;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    MultiplaneBuilder builder(cosm);
    std::ifstream lensdata(
        "tests/testdata.txt");

    if (!lensdata.is_open())
        throw(1);

    CompositeLensBuilder lensbuilder(z_d1);
    for (int i = 0; i < 9; i++) {
        float x, y;
        lensdata >> x >> y;
        x /= unit;
        y /= unit;
        Plummer plum(Dd1, 1e13 * MASS_SOLAR, 5, 1 / ANGLE_ARCSEC,
                     float2{.x = (float)x, .y = (float)y});
        lensbuilder.addLens(plum);
    }

    CompositeLensBuilder lensbuilder2(z_d2);
    for (int i = 0; i < 9; i++) {
        float x, y;
        lensdata >> x >> y;
        x /= unit;
        y /= unit;
        Plummer plum(Dd2, 1e13 * MASS_SOLAR, 6, 1 / ANGLE_ARCSEC,
                     float2{.x = (float)x, .y = (float)y});
        lensbuilder2.addLens(plum);
    }

    lensdata.close();

    builder.addPlane(lensbuilder);
    builder.addPlane(lensbuilder2);

    std::vector<float> srcs{0.2, 0.8, 1.5, 2.5};
    builder.setRedshifts(srcs);

    auto mp = builder.getCuMultiPlane();

    std::vector<std::vector<float>> factors{{0, 1, 0, 1, 0, 2, 0, 0.5, 0},
                                            {1, 1, 0, 0.5, 1.5, 1, 0, 0, 0.5}};
    mp.updateMassesCu(factors);

    float2 point{.x = -6.575343, .y = -6.692759};
    thrust::host_vector<float2> thetas;
    thetas.push_back(point);
    thrust::host_vector<float2> betas(1);

    thrust::device_vector<float2> d_thetas(thetas);
    thrust::device_vector<float2> d_betas(1);

    float2 *thetasptr = thrust::raw_pointer_cast(&d_thetas[0]);
    float2 *betasptr = thrust::raw_pointer_cast(&d_betas[0]);

    mp.traceThetas(thetasptr, betasptr, 1, 0);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x + 6.575343), 1e-5);
    EXPECT_LT(fabs(betas[0].y + 6.692759), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 1);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x + 0.071250), 1e-5);
    EXPECT_LT(fabs(betas[0].y + 1.502867), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 2);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x - 2.191279), 1e-5);
    EXPECT_LT(fabs(betas[0].y - 0.046877), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(thetasptr, betasptr, 1, 3);
    betas = d_betas;
    EXPECT_LT(fabs(betas[0].x - 2.964003), 1e-5);
    EXPECT_LT(fabs(betas[0].y - 0.463361), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.destroy();
}