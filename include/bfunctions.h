#pragma once
#include <vector>
#include "libslater.h"
#include <complex>
#include <iostream> //REMOVE - only added for my primitive caveman style debugging
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/legendre.hpp>
namespace slater {



class B_function_details {

private:
    double alpha; /// alpha of the radial part.
    Quantum_Numbers quantum_numbers; /// Set of quantum numbers for this function
    center_t center; /// the center of the function

public:
    /// Initialized details about one B functions
    /// \param exponent_ the exponnt used
    /// \param quantum_numbers
    /// \param center
    B_function_details(sto_exponent_t exponent_, const Quantum_Numbers &quantum_numbers,
        const center_t& center);

    sto_exponent_t get_alpha() const;
    const Quantum_Numbers &get_quantum_numbers() const;
    center_t get_center() const { return center; }

};


class B_functions_representation_of_STO {

    double b_functions_sum_rescaling = 0;
    std::vector<std::pair<double,B_function_details>> components;

    double calculate_coefficient(const Quantum_Numbers quantum_numbers,const unsigned int p) const;
public:
    B_functions_representation_of_STO(const STO_Basis_Function &sto, const center_t& new_center);
    auto size() { return components.size(); }
    auto begin() const { return components.cbegin(); }
    auto end() const { return components.cend(); }

    auto get_rescaling_coefficient() const { return b_functions_sum_rescaling; }
};


class B_function_Engine {
public:

    /// Return the value of a B-function with the given parameters at the point "r". A B-functions is completly
    /// defined by the three quantum numbers and alpha.
    /// \param quantum_numbers set of Quantum numbers, paramter
    /// \param alpha alpha parameter of the function
    /// \param r Point to calculate the value for
    /// \return value of the function at "r"
    std::complex<double> calculate(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &r) const;

    /// Evaluate Spherical Harmonics Y_l^m (theta,phi)
    std::complex<double> eval_spherical_harmonics(const Quantum_Numbers quantumNumbers,const double theta,const double phi) const;

    /// Converts cartesian representation of rr = (X,Y,Z) to spherical rr = (r,theta,phi)
    std::vector<double> cartesian_to_spherical(const center_t &r) const;
private:
};


}
