#include "sources.h"

#include <iostream>

uint8_t SourcePlane::check_hit(const Vector2D<float> &theta) const {
    // Check for each point |diff| < radius
	/*
    for (int i = 0; i < m_points_length; i++) {
        auto diff =
            sqrt((theta - m_points[i].position).lengthSq());
#ifndef __CUDA_ARCH__
        // std::cout << "Theta: " << theta.x() << "," << theta.y() << std::endl;
        // std::cout << "Diff: " << diff << std::endl;
#endif
        diff = fabs(diff);
#ifndef __CUDA_ARCH__
        // std::cout << "Diff: " << diff << std::endl;
#endif
        if (diff < m_points[i].radius) {
            return m_points[i].color;
        }
    }
	*/
    return 0;
}
