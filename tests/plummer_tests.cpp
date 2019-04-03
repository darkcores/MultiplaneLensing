#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

TEST(PlummerTests, TestBeta) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 1e13 * MASS_SOLAR, 30, 1 / ANGLE_ARCSEC,
                 Vector2D<float>(0, 0));
    std::vector<Vector2D<float>> betas;
    for (int i = 0; i < 10; i++) {
        Vector2D<float> point(i, i);
        auto trace = point - (plum.getAlpha(point) * (Dds / Ds));
        betas.push_back(trace);
    }

    EXPECT_LT(fabs(betas[0].x() - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[0].y() - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[1].x() - 0.9918526879), 1.25e-6);
    EXPECT_LT(fabs(betas[1].y() - 0.9918526879), 1.25e-6);
    EXPECT_LT(fabs(betas[2].x() - 1.9838130496), 1.25e-6);
    EXPECT_LT(fabs(betas[2].y() - 1.9838130496), 1.25e-6);
    EXPECT_LT(fabs(betas[3].x() - 2.9759840670), 1.25e-6);
    EXPECT_LT(fabs(betas[3].y() - 2.9759840670), 1.25e-6);
    EXPECT_LT(fabs(betas[4].x() - 3.9684597618), 1.25e-6);
    EXPECT_LT(fabs(betas[4].y() - 3.9684597618), 1.25e-6);
    EXPECT_LT(fabs(betas[5].x() - 4.9613217080), 1.25e-6);
    EXPECT_LT(fabs(betas[5].y() - 4.9613217080), 1.25e-6);
    EXPECT_LT(fabs(betas[6].x() - 5.9546365711), 1.25e-6);
    EXPECT_LT(fabs(betas[6].y() - 5.9546365711), 1.25e-6);
    EXPECT_LT(fabs(betas[7].x() - 6.9484547811), 1.25e-6);
    EXPECT_LT(fabs(betas[7].y() - 6.9484547811), 1.25e-6);
    EXPECT_LT(fabs(betas[8].x() - 7.9428103075), 1.25e-6);
    EXPECT_LT(fabs(betas[8].y() - 7.9428103075), 1.25e-6);
    EXPECT_LT(fabs(betas[9].x() - 8.9377213942), 1.25e-6);
    EXPECT_LT(fabs(betas[9].y() - 8.9377213942), 1.25e-6);
}

TEST(PlummerTests, TestMassUpdate) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 1e13 * MASS_SOLAR, 30, 1 / ANGLE_ARCSEC,
                 Vector2D<float>(0, 0));
    std::vector<Vector2D<float>> betas;
    plum.update(2.3);
    for (int i = 0; i < 10; i++) {
        Vector2D<float> point(i, i);
        auto trace = point - (plum.getAlpha(point) * (Dds / Ds));
        betas.push_back(trace);
    }
    EXPECT_LT(fabs(betas[0].x() - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[0].y() - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[1].x() - 0.9812611822), 1.25e-6);
    EXPECT_LT(fabs(betas[1].y() - 0.9812611822), 1.25e-6);
    EXPECT_LT(fabs(betas[2].x() - 1.9627700140), 1.25e-6);
    EXPECT_LT(fabs(betas[2].y() - 1.9627700140), 1.25e-6);
    EXPECT_LT(fabs(betas[3].x() - 2.9447633542), 1.25e-6);
    EXPECT_LT(fabs(betas[3].y() - 2.9447633542), 1.25e-6);
    EXPECT_LT(fabs(betas[4].x() - 3.9274574523), 1.25e-6);
    EXPECT_LT(fabs(betas[4].y() - 3.9274574523), 1.25e-6);
    EXPECT_LT(fabs(betas[5].x() - 4.9110399283), 1.25e-6);
    EXPECT_LT(fabs(betas[5].y() - 4.9110399283), 1.25e-6);
    EXPECT_LT(fabs(betas[6].x() - 5.8956641134), 1.25e-6);
    EXPECT_LT(fabs(betas[6].y() - 5.8956641134), 1.25e-6);
    EXPECT_LT(fabs(betas[7].x() - 6.8814459966), 1.25e-6);
    EXPECT_LT(fabs(betas[7].y() - 6.8814459966), 1.25e-6);
    EXPECT_LT(fabs(betas[8].x() - 7.8684637072), 1.25e-6);
    EXPECT_LT(fabs(betas[8].y() - 7.8684637072), 1.25e-6);
    EXPECT_LT(fabs(betas[9].x() - 8.8567592066), 1.25e-6);
    EXPECT_LT(fabs(betas[9].y() - 8.8567592066), 1.25e-6);
}
