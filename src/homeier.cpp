#include "homeier.h"
#include "logger.h"

namespace slater {

Homeier_Integrator::~Homeier_Integrator()
{

}

void Homeier_Integrator::init(const STO_Integration_Options &params)
{
    params.get(Use_Normalized_B_Functions_Parameter_Name, use_normalized_b_functions);
}

integral_value Homeier_Integrator::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    return 0;
}



}