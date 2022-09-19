#pragma once

#include "libslater.h"

namespace slater {

class Homeier_Integrator : public STO_Integrator {

    bool use_normalized_b_functions = 0;

public:



    virtual ~Homeier_Integrator();

    virtual void init(const STO_Integration_Options &params) ;
    virtual integral_value overlap(const std::array<STO_Basis_Function, 4> &) ;

};

}