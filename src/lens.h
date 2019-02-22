#pragma once

#include "util/vector2d.h"

template <class T> class Lens {
  private:
    Lens();

  public:
    virtual Vector2D<T> traceTheta(T Ds, T Dds, Vector2D<T> theta) const = 0;
    virtual Vector2D<T> getAlphaVector(Vector2D<T> theta) const = 0;
    virtual void getAlphaVectorDerivatives(Vector2D<T> theta, T &axx, T &ayy,
                                           T &axy) const = 0;
	virtual T getSurfaceMassDensity(Vector2D<T> theta) const = 0;
};
