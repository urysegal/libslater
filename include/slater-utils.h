#pragma once

#include "libslater.h"
#include "coordinates.h"

namespace slater {

/// Return the distance between two 3D points
/// \param A point A
/// \param B point B
/// \return distance between the points.
    double distance(
            const center_t &A,
            const center_t &B
    );

/// Return the magnitude of a vector
/// \param A vector
/// \return magnitude of vector.
    double vector_length(
            const center_t &A
    );

/// Return the vector connecting A to B so that A+Return=B
/// \param A first point
/// \param B second point
/// \return point R so that A+R=B
    center_t vector_between(
            const center_t &A,
            const center_t &B
    );


    center_t scale_vector(
            const center_t &P,
            double factor
    );

    std::complex<double> eval_spherical_harmonics(
            const Quantum_Numbers &quantumNumbers,
            const Spherical_Coordinates &spherical
    );
    std::complex<double> eval_spherical_harmonics_real(
            const Quantum_Numbers &quantumNumbers,
            const Spherical_Coordinates &spherical
    );

    double compute_reduced_bessel_function_half(
            const double order,
            const double z
    );
    double pochhammer(
            const double x,
            const double n
    );
}