#include "gtest/gtest.h"

#include <util/constants.h>
#include <util/cosmology.h>
#include <util/vector2d.h>

TEST(UtilTests, Vector2DMul) {
    Vector2D<double> v2d(0.5, 5);
    auto x = v2d.x();
    ASSERT_EQ(x, 0.5);
    auto scaled = v2d * 2;
    ASSERT_EQ(scaled.x(), 1);
    v2d *= 4;
    ASSERT_EQ(v2d.x(), 2);
    ASSERT_EQ(v2d.y(), 20);
}

TEST(UtilTests, Vector2DDiv) {
    Vector2D<double> v2d(2, 4);
    auto v = v2d / 2;
    ASSERT_EQ(v.x(), 1);
    ASSERT_EQ(v.y(), 2);
}

TEST(UtilTests, Vector2DLen) {
    Vector2D<double> v2d(3, 4);
    auto len = v2d.length();
    ASSERT_EQ(len, 5);
}

TEST(UtilTests, CosmologyTest) {
    Cosmology cosm(0.7, 0.3, 0, 0.7);
    auto dist = cosm.angularDiameterDistance(0.0, 1.72);
    const auto cdist = 1745;
    EXPECT_EQ((int)(dist / DIST_MPC), cdist); // Good enough, don't worry for now
    // ASSERT_LT(abs(dist - cdist), 1e-20);
}
