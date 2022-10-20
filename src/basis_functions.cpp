#include "libslater.h"
#include "coordinates.h"
#include "slater-utils.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/legendre.hpp>

namespace bm = boost::math;

namespace slater {


STO_Basis_Function_Info::STO_Basis_Function_Info( sto_exponent_t exponent_,
                                                 const Quantum_Numbers &quantum_numbers_) : exponent(exponent_),
        quantum_numbers(quantum_numbers_)
{
    auto alpha = exponent;
    auto n = quantum_numbers.n;
    normalization_coefficient = pow(pow(2.0*alpha,2*n+1),bm::factorial<double>(2*n));
    normalization_coefficient = sqrt(normalization_coefficient);
}

const Quantum_Numbers &STO_Basis_Function_Info::get_quantum_numbers() const
{
    return quantum_numbers;
}

sto_exponent_t STO_Basis_Function_Info::get_exponent() const
{
    return exponent;
}

void STO_Basis_Function_Info::set_exponent(sto_exponent_t e)
{
    exponent = e;
}

sto_coefficient_t STO_Basis_Function_Info::get_coefficient() const
{
    return normalization_coefficient;
}


void STO_Basis_Function_Info::set_coefficient(sto_coefficient_t c)
{
    normalization_coefficient = c;
}


void STO_Basis_Function_Info::set_quantum_numbers(const Quantum_Numbers &quantumNumbers)
{
    quantum_numbers = quantumNumbers;
}

STO_Basis_Function::STO_Basis_Function(STO_Basis_Function_Info function_info_, center_t location_):
    function_info(function_info_), center(location_)

{}


const Quantum_Numbers &STO_Basis_Function::get_quantum_numbers() const
{
    return function_info.get_quantum_numbers();
}

sto_exponent_t STO_Basis_Function::get_exponent() const {
    return function_info.get_exponent();
}

sto_coefficient_t STO_Basis_Function::get_coefficient() const
{
    return function_info.get_coefficient();
}

center_t STO_Basis_Function::get_center() const
{
    return center;
}

std::complex<double> STO_Basis_Function::evaluate(const center_t &point) const
{
    // STO_{n,l}^m(alpha,rr) = prefactor * r^(n-1)*exp(-apha*r) * Y
    // rr = (r,theta,phi)
    // prefactor = sqrt((2*alpha)^(2n+1)/((2*n)!)
    // Y = Y_{l,m}(theta,phi)

    //COMPLEX ARITHMETIC NEEDS TO BE Checked
    auto n = this->get_quantum_numbers().n;
    auto alpha = this->get_exponent();

    // Cartesian Representation of r to Spherical representation
    Spherical_Coordinates spherical(point);

    //need to extract theta and phi from r_spherical
    auto Y = eval_spherical_harmonics(this->get_quantum_numbers(),spherical.theta,spherical.phi);

    return this->get_coefficient() * pow(spherical.radius,n-1)*exp(-alpha*spherical.radius) * Y;
}

std::complex<double> STO_Basis_Function::evaluate_conjugate(const center_t &point) const
{
    return std::conj(evaluate(point));
}


}
