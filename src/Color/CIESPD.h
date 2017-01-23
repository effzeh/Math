#pragma once

#define LIB_MATH_COLOR_CIE_SPD_H 1

#include "../Vector.h"

namespace math {

// Mean relative spectral radiant power distribution S0 and first two
// eigenvectors S1 and S2.
//
// Used in the CIE method of calculating daylight illuminants.
class CIESPD
{
    static const Vec3 Table[];

public:
    static Vec3 eval(float lambda)
    {
        if (lambda < 290 || lambda >= 840)
            return {0,0,0};

        auto i = static_cast<int>(Floor(lambda)) / 10;
        auto s = (lambda - 10 * i) / 10.0f;

        i -= 290 / 10;

        return Lerp(Table[i], Table[i + 1], s);
    }

    Vec3 operator()(float lambda) const
    {
        return eval(lambda);
    }
};

} // namespace math
