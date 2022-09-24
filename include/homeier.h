#pragma once

#include "libslater.h"
#include "bfunctions.h"

namespace slater {

class Homeier_Integrator : public STO_Integrator {

private:
    bool use_normalized_b_functions = 0;
    int number_of_quadrature_points = 1024;
    std::vector<std::pair<B_function_details, B_function_details> > equivalence_series;

    double calculate_one_weight(const B_function_details &f1, const B_function_details &f2);
    energy_unit_t integrate_with_b_functions(const B_function_details &f1, const B_function_details &f2) const;

    void create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2) ;
    double calculate_guassian_point(const B_function_details &f1, const B_function_details &f2, double s) const;

    double calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s) const ;
    double calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const;


public:



    virtual ~Homeier_Integrator();

    virtual void init(const STO_Integration_Options &params) override;

    virtual energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &) override;

};

}