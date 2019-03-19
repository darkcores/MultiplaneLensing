#include "sources.h"
#include "util/error.h"

SourcePlane SourcePlaneBuilder::getPlane() {
    size_t size = sizeof(SourceData) * m_points.size();
    if (size == 0) {
        std::cerr << "No sources added" << std::endl;
        std::terminate();
    }
    ptr = (SourceData *)malloc(size);
	cpuErrchk(ptr);
    memcpy(ptr, &m_points[0], size);

    return SourcePlane(m_redshift, ptr, m_points.size());
}
