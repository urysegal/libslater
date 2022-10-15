#pragma once

#include "libslater.h"

namespace slater {

/// Return the distance between two 3D points
/// \param A point A
/// \param B point B
/// \return distance between the points.
double distance(const center_t &A, const center_t &B) ;

}