#pragma once

#include "util/vector2d.h"

class Lens {
  private:

  public:
    virtual Vector2D<double> traceTheta(double Ds, double Dds,
                                             Vector2D<double> theta) const = 0;
    virtual Vector2D<double> getAlphaVector(Vector2D<double> theta) const = 0;
    virtual void getAlphaVectorDerivatives(Vector2D<double> theta, double &axx,
                                           double &ayy, double &axy) const = 0;
    virtual double getSurfaceMassDensity(Vector2D<double> theta) const = 0;
};
