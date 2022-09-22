#include "libslater.h"

namespace slater {


STO_Basis_Function_Info::STO_Basis_Function_Info(sto_coefficient_t coefficient_, sto_exponent_t exponent_,
                                                 moment_quantum_number_t moment_quantum_number_,
                                                 principal_quantum_number_t principal_quantum_number_,
                                                 orientation_quantum_number_t orientation_) :
                                                 coefficient(coefficient_), exponent(exponent_),
                                                 moment_quantum_number(moment_quantum_number_),
                                                 principal_quantum_number(principal_quantum_number_),
                                                 orientation(orientation_)
{}

STO_Basis_Function::STO_Basis_Function(STO_Basis_Function_Info function_info_, center_t location_):
    function_info(function_info_), center(location_)

{}

}
