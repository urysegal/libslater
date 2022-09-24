#include "bfunctions.h"

namespace slater
{


B_function_details::B_function_details(sto_coefficient_t coefficient_, sto_exponent_t exponent_,
                                       const Quantum_Numbers &quantum_numbers_, const center_t &center_) :
                                       coefficient(coefficient_), exponent(exponent_),
                                       quantum_numbers(quantum_numbers_), center(center_)
{}


sto_coefficient_t B_function_details::get_coefficient() const
{
    return coefficient;
}

sto_exponent_t B_function_details::get_exponent() const
{
    return exponent;
}

const Quantum_Numbers &B_function_details::get_quantum_numbers() const
{
    return quantum_numbers;
}



B_functions_representation_of_STO::B_functions_representation_of_STO(const STO_Basis_Function &sto)
{
    /// Gautam: create all the necessary B functions from the given sto. Here I just make up two.

    {
        B_function_details bfunc(sto.get_coefficient()/7, sto.get_exponent()*0.98, sto.get_quantum_numbers(), sto.get_center());
        components.emplace_back(bfunc);
    }

    {
        // Gautam - You can play with the quantum numbers too.
        auto quantum_numbers = sto.get_quantum_numbers();
        quantum_numbers.l++;

        B_function_details bfunc(sto.get_coefficient()*0.9, sto.get_exponent()/2, quantum_numbers, sto.get_center());
        components.emplace_back(bfunc);
    }
}

}
