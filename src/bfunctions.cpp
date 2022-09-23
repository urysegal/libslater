#include "bfunctions.h"

namespace slater
{

B_functions_representation_of_STO::B_functions_representation_of_STO(const STO_Basis_Function &sto)
{
        /// Gotam: create all the necessary B functions from the given sto.

        {
            B_function_details bfunc(sto.get_coefficient(), sto.get_exponent(), sto.get_quantum_numbers(), sto.get_center());
            components.emplace_back(bfunc);
            this->coefficient.emplace_back(calculate_coefficient(bfunc));
        }

    {
        // Gotham - You can play with the quantum numbers too.
        auto quantum_numbers = sto.get_quantum_numbers();
        quantum_numbers.l++;

        B_function_details bfunc(sto.get_coefficient()*0.9, sto.get_exponent()/2, quantum_numbers, sto.get_center());
        components.emplace_back(bfunc);
        this->coefficient.emplace_back(calculate_coefficient(bfunc));
    }


}

}
