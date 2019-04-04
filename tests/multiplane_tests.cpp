#include "gtest/gtest.h"

#include <multiplane.h>

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
                         Vector2D<float>(x, y));
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

    auto mp = builder.getMultiPlane();

    Vector2D<float> point(2, 1);
    std::vector<Vector2D<float>> thetas{point};
    std::vector<Vector2D<float>> betas(1);

    mp.traceThetas(&thetas[0], &betas[0], 1, 0);
    EXPECT_EQ(betas[0].x(), 2.0);
    EXPECT_EQ(betas[0].y(), 1.0);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 1);
    EXPECT_LT(fabs(betas[0].x() + 1.579611), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 0.840784), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 2);
    EXPECT_LT(fabs(betas[0].x() + 1.736920), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 0.876505), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 3);
    EXPECT_LT(fabs(betas[0].x() + 1.310565), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 0.621893), 1e-5);
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

    auto mp = builder.getMultiPlane();

    Vector2D<float> point(2, 1);
    std::vector<Vector2D<float>> thetas{point};
    std::vector<Vector2D<float>> betas(1);

    std::vector<std::vector<float>> factors{{0, 1, 0, 1, 0, 0.5, 0, 2, 0},
                                            {1, 0.5, 0, 1, 1.5, 0, 0, 1, 0.5}};
    mp.updateMasses(factors);

    mp.traceThetas(&thetas[0], &betas[0], 1, 0);
    EXPECT_EQ(betas[0].x(), 2.0);
    EXPECT_EQ(betas[0].y(), 1.0);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 1);
    EXPECT_LT(fabs(betas[0].x() - 4.607982), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 1.687718), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 2);
    EXPECT_LT(fabs(betas[0].x() - 3.180912), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 2.014467), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.traceThetas(&thetas[0], &betas[0], 1, 3);
    EXPECT_LT(fabs(betas[0].x() - 1.663422), 1e-5);
    EXPECT_LT(fabs(betas[0].y() + 1.857666), 1e-5);
    // printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.destroy();
}
