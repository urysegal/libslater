#pragma once

#include "libslater.h"
#include "bfunctions.h"

namespace slater {

class Homeier_Integrator : public STO_Integrator {

private:
    bool use_normalized_b_functions = 0;
    int number_of_quadrature_points = 1024;
    std::vector<std::pair<B_function_details&, B_function_details&> > components;

    double calculate_one_weight(const B_function_details &f1, const B_function_details &f2);


public:



    virtual ~Homeier_Integrator();

    virtual void init(const STO_Integration_Options &params) ;

    virtual integral_value overlap(const std::array<STO_Basis_Function, 2> &) override;

};

}