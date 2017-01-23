#pragma once

#define LIB_MATH_COLOR_CIE_DSERIES_H 1

#include "CIESPD.h"

namespace math {

// Computes the relative spectral power distribution (SPD) of a D series
// illuminant.
// Normalized such that P(560) = 1.
class CIEDSeries
{
    // https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

    float M1;
    float M2;

public:
    CIEDSeries() : CIEDSeries(6504.0f)
    {
    }

    explicit CIEDSeries(float t)
    {
        t = Clamp(t, 4000.0f, 25000.0f);

        float x = (t <= 7000.0f)
            ? (0.244063f + 0.09911e+3f / t + 2.9678e+6f / (t*t) - 4.6070e+9f / (t*t*t))
            : (0.237040f + 0.24748e+3f / t + 1.9018e+6f / (t*t) - 2.0064e+9f / (t*t*t));
        float y = -3.0f * (x*x) + 2.87f * x - 0.275f;
        float M = 0.0241f + 0.2562f * x - 0.7341f * y;

        M1 = (-1.3515f -  1.7703f * x +  5.9114f * y) / M;
        M2 = ( 0.0300f - 31.4424f * x + 30.0717f * y) / M;
    }

    // Evaluate at the given wavelength.
    // Lambda is in nm.
    float operator()(float lambda) const
    {
        auto s = CIESPD::eval(lambda);
        return s.x + M1 * s.y + M2 * s.z;
    }
};

} // namespace math
