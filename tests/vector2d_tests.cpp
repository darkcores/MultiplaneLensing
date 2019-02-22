#include "gtest/gtest.h"

#include <util/vector2d.h>

TEST(Vector2DTests, ScalingDouble) {
	Vector2D<double> v2d(0.5, 5);
	auto x = v2d.x();
	ASSERT_EQ(x, 0.5);
	auto scaled = v2d * 2;
	ASSERT_EQ(scaled.x(), 1);
	v2d *= 4;
	ASSERT_EQ(v2d.x(), 2);
	ASSERT_EQ(v2d.y(), 20);
}
