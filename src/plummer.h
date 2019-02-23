#pragma once

#include "lens.h"
// #include "util/vector2d.h"

class Plummer : public Lens {
  private:
	double m_Dd;
    double m_mass;
    double m_angularwidth, m_angularwidth2;
    Vector2D<double> m_angularpos;

  public:
    Plummer(const double Dd, const double mass, const double angularwidth,
            const Vector2D<double> &angularposition);

    virtual Vector2D<double> traceTheta(double Ds, double Dds,
                                Vector2D<double> theta) const;
    virtual Vector2D<double> getAlphaVector(Vector2D<double> theta) const;
    virtual void getAlphaVectorDerivatives(Vector2D<double> theta, double &axx,
                                   double &ayy, double &axy) const;
    virtual double getSurfaceMassDensity(Vector2D<double> theta) const;


  private:
    double getMassInside(double thetaLength) const;
	double getLensDistance() const	{ return m_Dd; }
	double getProfileSurfaceMassDensity(double thetaLength) const;
};
