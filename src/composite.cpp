#include "composite.h"
#include "util/error.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

void CompositeLensBuilder::addLens(const Plummer &lens) {
    m_lenses.push_back(lens);
}

CompositeLens CompositeLensBuilder::getLens() {
    // m_lenses.push_back(LensData());
	Plummer *lens_ptr = nullptr;
    size_t size = sizeof(Plummer) * m_lenses.size();
    if (size == 0) {
        std::cerr << "No lenses added" << std::endl;
        std::terminate();
    }
    lens_ptr = (Plummer *)malloc(size);
    cpuErrchk(lens_ptr);
    std::memcpy((void *)lens_ptr, &m_lenses[0], size);
    return CompositeLens(lens_ptr, m_lenses.size());
}
