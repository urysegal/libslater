#include "libslater.h"

namespace slater {


STO_Basis_Function_Info::STO_Basis_Function_Info(sto_coefficient_t coefficient_, sto_exponent_t exponent_,
                                                 const Quantum_Numbers &quantum_numbers_) :
                                                 coefficient(coefficient_), exponent(exponent_),
                                                 quantum_numbers(quantum_numbers_)
{}

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
    return coefficient;
}


void STO_Basis_Function_Info::set_coefficient(sto_coefficient_t c)
{
    coefficient = c;
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

}
