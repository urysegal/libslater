#include "libslater.h"
#include "coordinates.h"
#include "slater-utils.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>

namespace bm = boost::math;

namespace slater {


void Quantum_Numbers::validate() const
{
    assert(l < n);
    assert(l >= 0);
    assert(abs(m) <= l );
}

STO_Basis_Function_Info::STO_Basis_Function_Info( sto_exponent_t exponent_,
                                                 const Quantum_Numbers &quantum_numbers_) : exponent(exponent_),
        quantum_numbers(quantum_numbers_)
{
    auto alpha = exponent;
    auto n = int(quantum_numbers.n);
    normalization_coefficient = pow(2.0*alpha,2*n+1)/bm::factorial<double>(2*n);
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

sto_coefficient_t STO_Basis_Function::get_normalization_coefficient() const
{
    return function_info.get_coefficient();
}

center_t STO_Basis_Function::get_center() const
{
    return center;
}

std::complex<double> STO_Basis_Function::evaluate(const center_t &point) const
{
    auto n = this->get_quantum_numbers().n;
    auto alpha = this->get_exponent();

    Spherical_Coordinates spherical(point);

    auto Y = eval_spherical_harmonics(this->get_quantum_numbers(),spherical.theta,spherical.phi);

    return this->get_normalization_coefficient() * pow(spherical.radius, n - 1) * exp(-alpha * spherical.radius) * Y;
}

std::complex<double> STO_Basis_Function::evaluate_conjugate(const center_t &point) const
{
    return std::conj(evaluate(point));
}


}
