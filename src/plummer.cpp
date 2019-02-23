#include "plummer.h"

#include "util/constants.h"
#include <cmath>

Plummer::Plummer(const double Dd, const double mass, const double angularwidth,
                 const Vector2D<double> &angularposition) {
    m_Dd = Dd;
    m_mass = mass;
    m_angularwidth = angularwidth;
    m_angularwidth2 = angularwidth * angularwidth;
    m_angularpos = angularposition;
}

Vector2D<double> Plummer::traceTheta(double Ds, double Dds,
                                     Vector2D<double> theta) const {
    auto trace = getAlphaVector(theta);
    trace = theta - (trace * (Dds / Ds));
    return trace;
}

Vector2D<double> Plummer::getAlphaVector(Vector2D<double> theta) const {
    Vector2D<double> alphavect(0, 0);
    double mass = getMassInside(theta.length());

    if (mass != 0) {
        double massScale = fabs(mass);
        double scaledMass = mass / massScale;
        double scalefact = sqrt((4.0 * CONST_G * massScale) /
                                (SPEED_C * SPEED_C * getLensDistance()));
        Vector2D<double> scaledtheta = theta / scalefact;

        double length2 = scaledtheta.lengthSq();

        if (length2 != 0) // no need to divide if theta == (0, 0)
            alphavect = (theta / (scaledtheta.lengthSq())) * scaledMass;
    }
    return alphavect;
}

double Plummer::getMassInside(double thetaLength) const {
    double t = thetaLength / m_angularwidth;
    double x = t * t;
    double x2 = 1.0 / x;

    return m_mass * (1.0 / (1.0 + x2));
}
void Plummer::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx,
                                        double &ayy, double &axy) const {
    double axx0 = 0;
    double axy0 = 0;
    double ayy0 = 0;
    double thetaLength = theta.length();
    double massInside = getMassInside(thetaLength);

    if (massInside != 0) {
        double massScale = fabs(massInside);
        double scaledMass = massInside / massScale;
        double scaleFactor = sqrt(4.0 * CONST_G * massScale /
                                  (SPEED_C * SPEED_C * getLensDistance()));
        Vector2D<double> scaledTheta = theta / scaleFactor;
        double t2 = scaledTheta.lengthSq();
        double t4 = t2 * t2;
        double x = scaledTheta.x();
        double y = scaledTheta.y();
        double x2 = x * x;
        double y2 = y * y;
        double xy = x * y;
        double part1 = 0;
        double part3 = 0;

        if (t4 != 0) {
            part1 = scaledMass * (y2 - x2) / t4;
            part3 = -scaledMass * 2.0 * xy / t4;
        }

        axx0 += part1;
        ayy0 += -part1;
        axy0 += part3;
    }

    double factor2 = (8.0 * CONST_PI * CONST_G * getLensDistance() *
                      getProfileSurfaceMassDensity(thetaLength)) /
                     (SPEED_C * SPEED_C);
    double thetaLengthSquared = theta.lengthSq();

    if (thetaLengthSquared != 0) {
        double x = theta.x();
        double y = theta.x();
        double x2 = x * x;
        double y2 = y * y;
        double xy = x * y;

        double part2 = factor2 / thetaLengthSquared;

        axx0 += part2 * x2;
        ayy0 += part2 * y2;
        axy0 += part2 * xy;
    } else {
        axx0 += factor2 * 0.5;
        ayy0 += factor2 * 0.5;
        // third part should be zero in (0,0) by symmetry (choice of axis
        // orientation doesn't change anything about the view)
    }

    axx = axx0;
    ayy = ayy0;
    axy = axy0;

    // return true;
}

double Plummer::getProfileSurfaceMassDensity(double thetaLength) const {
    double scaledtheta = thetaLength / m_angularwidth;
    double denom = (1.0 + scaledtheta * scaledtheta) * getLensDistance();
    double dens = m_mass / (CONST_PI * m_angularwidth2 * denom * denom);
    return dens;
}

double Plummer::getSurfaceMassDensity(Vector2D<double> theta) const {
    return getProfileSurfaceMassDensity(theta.length());
}
