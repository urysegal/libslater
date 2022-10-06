#pragma once
#include "libslater.h"

namespace slater {

struct Spherical_Coordinates {
    double theta;
    double phi;
    double radius;
    Spherical_Coordinates(const center_t &cartesian);
};


}