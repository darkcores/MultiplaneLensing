#pragma once

#include <cstdio>
#include <cstdlib>

#ifdef __CUDACC__

#define gpuErrchk(ans)                                                         \
    { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
                line);
        if (abort)
            throw(code);
    }
}

#endif

#define cpuErrchk(ans)                                                         \
    { cpuAssert((ans), __FILE__, __LINE__); }
inline void cpuAssert(void *code, const char *file, int line,
                      bool abort = true) {
    if (code == NULL) {
        fprintf(stderr, "Assert memory error: %s %d\n", file, line);
        if (abort)
            throw(-1);
    }
}
