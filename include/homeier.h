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

    /// Build an Homeier-style integrator
    Homeier_Integrator();

    /// Release any memory used by the integrator
    virtual ~Homeier_Integrator();

    /// Initialize the Homeier integrator with a set of options
    /// \param params set of options for the integrator
    virtual void init(const STO_Integration_Options &params) override;

    /// Calculate the Overlap integral between the two given STO basis functions
    /// \return the resulting energy quantity
    virtual energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &) override;


private:

    bool use_normalized_b_functions = 0; /// Should we calculate with normalized B functions?
    int number_of_quadrature_points = 1024; /// How many quadrature points we should calculate

    /// When integrating of each pair of functions in this vector and then summing up the values, you get the
    /// Overlap integral value
    std::vector<std::pair<
    std::pair<double,B_function_details>,
    std::pair<double, B_function_details> > > equivalence_series;

    /// Create all the pair of B functions and their normalization coefficients from the two given sequences of B functions,
    /// each representing an STO. The result is kept in the "equivalence_series" member
    /// \param f1 First sequences of B functions
    /// \param f2 Second sequences of B functions
    void create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2) ;

    /// Fiven two specific B functions, calculate the overlap integral
    /// \param f1 First B function
    /// \param f2 Second B function
    /// \return Partial overlap integral value
    double integrate_overlap_using_b_functions(const B_function_details &f1, const B_function_details &f2) const;

    /// This function is called back from the Gaussian Quadrature mechanism to get one value, at s, of the overlap
    /// integral.
    /// \param f1 details of the first B functions
    /// \param f2 details of the second B functions
    /// \param s point at which to calculate the integral
    /// \return overlap integral value at s
    double calculate_overlap_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const;



    double calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s) const ;
    double calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const;

};

}