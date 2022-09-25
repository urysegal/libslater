#pragma once
#include <vector>
#include "libslater.h"

namespace slater {



class B_function_details {

private:
    sto_coefficient_t coefficient; /// Normalization Coefficient of the B function
    sto_exponent_t exponent; /// exponent of the radial part.
    Quantum_Numbers quantum_numbers; /// Set of quantum numbers for this function
    center_t center; /// the center of the function

public:
    /// Initialized details about one B functions
    /// \param coefficient_ a coefficient for the function, if desired. Use 1 if irrelevant for your case
    /// \param exponent_ the exponnt used
    /// \param quantum_numbers
    /// \param center
    B_function_details(sto_coefficient_t coefficient_, sto_exponent_t exponent_, const Quantum_Numbers &quantum_numbers,
        const center_t& center);

    sto_coefficient_t get_coefficient() const ;
    sto_exponent_t get_exponent() const;
    const Quantum_Numbers &get_quantum_numbers() const;


};


class B_functions_representation_of_STO {

    std::vector<B_function_details> components;

    double calculate_coefficient(const B_function_details &) const;
public:
    B_functions_representation_of_STO(const STO_Basis_Function &sto);
    auto size() { return components.size(); }
    auto begin() const { return components.cbegin(); }
    auto end() const { return components.cend(); }
};




}
