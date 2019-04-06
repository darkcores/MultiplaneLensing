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
	// printf("Composite size: %lu\n", sizeof(CompositeLens));
	// printf("Plummer size: %lu\n", sizeof(Plummer)); 
    Plummer *lens_ptr;
    size_t size = m_lenses.size();
	// printf("Lenses in composite: %lu\n", size);
    if (size == 0) {
        std::cerr << "No lenses added" << std::endl;
        throw(-1);
    }
    size_t numbytes = sizeof(Plummer) * size;
	// printf("Bytes: %lu\n", numbytes);
    gpuErrchk(cudaMalloc(&lens_ptr, numbytes));
    gpuErrchk(
        cudaMemcpy(lens_ptr, &m_lenses[0], numbytes, cudaMemcpyHostToDevice));
    return CompositeLens(lens_ptr, (int)size, true);
}
