#include "gtest/gtest.h"

#include <util/vector2d.h>
#include <util/cosmology.h>

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

TEST(UtilTests, CosmologyTest) {
	Cosmology cosm(0.7, 0.3, 0, 0.7);
	auto dist = cosm.angularDiameterDistance(0.0, 1.72);
	ASSERT_EQ((int)dist, 1745); // Good enough, don't worry for now
}
