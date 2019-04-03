#include "composite.h"
#include "util/error.h"
#include <iostream>

int CompositeLens::destroy() {
    if (m_cuda) {
        gpuErrchk(cudaFree(m_lenses));
    } else {
        free(m_lenses);
    }
    return 0;
}

CompositeLens CompositeLensBuilder::getCuLens() {
    Plummer *lens_ptr;
#ifdef __CUDACC__
    size_t size = m_lenses.size();
    if (size == 0) {
        std::cerr << "No lenses added" << std::endl;
        throw(-1);
    }
    size_t numbytes = sizeof(Plummer) * size;
    gpuErrchk(cudaMalloc(&lens_ptr, numbytes));
    gpuErrchk(
        cudaMemcpy(lens_ptr, &m_lenses[0], numbytes, cudaMemcpyHostToDevice));
#else
    CompositeLens *lens_ptr = nullptr;
#endif
    return CompositeLens(lens_ptr, size, true);
}
