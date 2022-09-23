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

void STO_Basis_Function_Info::set_quantum_numbers(const Quantum_Numbers &quantumNumbers)
{
    quantum_numbers = quantumNumbers;
}

STO_Basis_Function::STO_Basis_Function(STO_Basis_Function_Info function_info_, center_t location_):
    function_info(function_info_), center(location_)

{}

}
