#pragma once

#include "libslater.h"
#include "bfunctions.h"

namespace slater {


/// This class implements the calculation of the Overlap integral between two STOs using B functions, using
/// this work:
/// On the evaluation of overlap integrals with exponential‐type basis functions, HHH Homeier, EO Steinborn,
/// International journal of quantum chemistry 42 (4), 761-778

class Homeier_Integrator : public STO_Integrator {

public:

    /// Release any memory used by the integrator
    virtual ~Homeier_Integrator();

    virtual void init(const STO_Integration_Options &params) override;
    virtual energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &) override;


private:

    bool use_normalized_b_functions = 0; /// Should we calculate with normalized B functions?
    int number_of_quadrature_points = 1024; /// How many quadrature points we should calculate

    /// When integrating of each pair of functions in this vector and then summing up the values, you get the
    /// Overlap integral value
    std::vector<std::pair<B_function_details, B_function_details> > equivalence_series;

    void create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2) ;

    energy_unit_t integrate_with_b_functions(const B_function_details &f1, const B_function_details &f2) const;

    double calculate_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const;

    double calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s) const ;
    double calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const;



};

}