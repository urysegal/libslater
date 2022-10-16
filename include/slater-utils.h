#pragma once

#include "libslater.h"

namespace slater {

/// Return the distance between two 3D points
/// \param A point A
/// \param B point B
/// \return distance between the points.
double distance(
    const center_t &A,
    const center_t &B
);


/// Return the vector connecting A to B so that A+Return=B
/// \param A first point
/// \param B second point
/// \return point R so that A+R=B
center_t vector_between(
    const center_t &A,
    const center_t &B
);


std::complex<double> eval_spherical_harmonics(
    const Quantum_Numbers quantumNumbers,
    const double theta,
    const double phi
);



}