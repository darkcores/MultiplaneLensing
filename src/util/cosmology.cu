#include "cosmology.h"

#include "constants.h"
#include <cmath>

__host__ __device__ Cosmology::Cosmology(const double h, const double omega_m,
                                         const double omega_r,
                                         const double omega_v, const double w) {
    m_h = h;
    m_Wm = omega_m;
    m_Wr = omega_r;
    m_Wv = omega_v;
    m_Wk = 1 - omega_m - omega_r - omega_v;
    m_w = w;
}

__host__ __device__ double
Cosmology::angularDiameterDistance(double z1, double z2, int num) const {
    if (z2 < 0.0) {
        z2 = z1;
        z1 = 0;
    }

    if (z1 > z2) {
        auto tmp = z1;
        z1 = z2;
        z2 = tmp;
    }

    if (z1 < 0 || z2 < 0) {
        // TODO actual error handling
        return -1.0;
    }

    double sum = 0;
    double R2 = 1.0 / (1.0 + z1);
    double R1 = 1.0 / (1.0 + z2);
    double dR = (R2 - R1) / ((double)num);
    double R = R1;
    for (int i = 0; i < num; i++, R += dR) {
        double term1 = m_Wv * pow(R, 1.0 - 3.0 * m_w);
        double term2 = m_Wm * R;
        double term3 = m_Wr;
        double term4 = m_Wk * R * R;

        double val1 = 1.0 / sqrt(term1 + term2 + term3 + term4);

        term1 = m_Wv * pow(R + dR, 1.0 - 3.0 * m_w);
        term2 = m_Wm * (R + dR);
        term3 = m_Wr;
        term4 = m_Wk * (R + dR) * (R + dR);

        double val2 = 1.0 / sqrt(term1 + term2 + term3 + term4);

        sum += ((val1 + val2) / 2.0) * dR; // trapezium rule
    }

    double A = 0;

    if (m_Wk == 0)
        A = sum;
    else if (m_Wk > 0)
        A = (1.0 / sqrt(m_Wk)) * sinh(sqrt(m_Wk) * sum);
    else // W_k < 0
        A = (1.0 / sqrt(-m_Wk)) * sin(sqrt(-m_Wk) * sum);

    double result = SPEED_C * A * ((1.0 / (1.0 + z2)) * TH) / m_h;

    return result;
}

__host__ __device__ double
Cosmology::redshiftForAngularDiameterDistance(double Dtarget, double zref,
                                              double zmax) const {
    return -1.0;
}
