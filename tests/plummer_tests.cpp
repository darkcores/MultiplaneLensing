#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

TEST(PlummerTests, TestBetaF) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.5;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 100000000000000, 30 * ANGLE_ARCSEC);
    plum.setSource(Ds, Dds);
    plum.setScale(60 * 60);
    Vector2D<float> point(1 * ANGLE_ARCMIN, 2 * ANGLE_ARCMIN);
    point *= (60 * 60);
    Vector2D<double> rpoint(1 * ANGLE_ARCMIN, 2 * ANGLE_ARCMIN);
    auto trace = plum.getBetaf(point);
    trace /= (60 * 60);
    // auto rtrace = plum.getBeta(rpoint);
    // EXPECT_EQ(trace.x(), rtrace.x());
    // EXPECT_EQ(trace.y(), rtrace.y());
    EXPECT_LT(abs(trace.x() - 0.00029089), 1e-8);
    EXPECT_LT(abs(trace.y() - 0.00058178), 1e-8);
}

// To look for a wierd result in the API
TEST(PlummerTests, TestBeta11) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 1e13 * MASS_SOLAR * 11, 30 * ANGLE_ARCSEC);
    plum.setSource(Ds, Dds);
    plum.setScale(60 * 60);
    Vector2D<float> point(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC);
    point *= (60 * 60);
    auto trace = plum.getBetaf(point);
    trace /= (60 * 60);
    // EXPECT_EQ(trace.x(), 4.41364469e-6);
    // EXPECT_EQ(trace.y(), 4.41364469e-6);
    EXPECT_LT(abs(trace.x() - 4.41364469e-6), 1e-12);
    EXPECT_LT(abs(trace.y() - 4.41364469e-6), 1e-12);
}

TEST(PlummerTests, TestMassUpdate) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 1e13 * 0, 30 * ANGLE_ARCSEC);
	plum.setMass(1e13 * MASS_SOLAR * 11);
    plum.setSource(Ds, Dds);
    plum.setScale(60 * 60);
    Vector2D<float> point(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC);
    point *= (60 * 60);
    auto trace = plum.getBetaf(point);
    trace /= (60 * 60);
    // EXPECT_EQ(trace.x(), 4.41364469e-6);
    // EXPECT_EQ(trace.y(), 4.41364469e-6);
    EXPECT_LT(abs(trace.x() - 4.41364469e-6), 1e-12);
    EXPECT_LT(abs(trace.y() - 4.41364469e-6), 1e-12);
}
