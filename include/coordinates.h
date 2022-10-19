#pragma once
#include "libslater.h"

namespace slater {

struct Spherical_Coordinates {
    double theta = 0;
    double phi = 0;
    double radius = 0;
    Spherical_Coordinates(const center_t &cartesian);
    Spherical_Coordinates() = default;
};


}