#pragma once
#include <vector>
#include "libslater.h"

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

    double calculate_coefficient(const B_function_details &bfd,const unsigned int p) const;
public:
    B_functions_representation_of_STO(const STO_Basis_Function &sto, const center_t& new_center);
    auto size() { return components.size(); }
    auto begin() const { return components.cbegin(); }
    auto end() const { return components.cend(); }

    auto get_rescaling_coefficient() const { return b_functions_sum_rescaling; }
};




}
