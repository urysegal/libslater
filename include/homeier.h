#pragma once

#include "libslater.h"
#include "bfunctions.h"
#include "gaunt.h"
#include <boost/math/special_functions/factorials.hpp>
#include "integrators.h"


namespace slater {

const std::string overlap_homeier_imp_name = "overlap-B-functions-homeier";
const std::string kinetic_homeier_imp_name = "kinetic-B-functions-homeier";


void shift_first_center_to_origin(const center_t &c1, center_t c2, center_t *new_centers);


/// This class implements the calculation of the Overlap integral between two STOs using B functions, using
/// this work:
/// On the evaluation of overlap integrals with exponential‐type basis functions, HHH Homeier, EO Steinborn,
/// International journal of quantum chemistry 42 (4), 761-778

class Homeier_Integrator : public STO_Integrator {

public:

    /// Build an Homeier-style integrator
    explicit Homeier_Integrator(integration_types t) :
            STO_Integrator(2,1), which_type(t)
    {}

    /// This constructor is only used for building the factory
    explicit Homeier_Integrator(const std::string &name) noexcept : STO_Integrator(name)
    {
        if ( name == kinetic_homeier_imp_name ) {
            which_type = integration_types::KINETIC;
        } else {
            which_type = integration_types::OVERLAP;
        }
    }

    [[nodiscard]] STO_Integrator *clone() const override;

    /// Initialize the Homeier integrator with a set of options
    /// \param params set of options for the integrator
    void init(const STO_Integration_Options &params) override;

    energy_unit_t integrate(
            const std::vector<STO_Basis_Function> &functions,
            const std::vector<center_t> &centers
    ) override;

    static double calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s)  ;
    static std::complex<double> calculate_S(const B_function_details &f1, const B_function_details &f2, double s) ;

    energy_unit_t calculate_B_function_kinetic(B_function_details p1,
                                               B_function_details p2,
                                               energy_unit_t &S1,
                                               energy_unit_t &S2,
                                               bool normalized);


private:

    integration_types which_type = integration_types::OVERLAP;

    static B_function_Engine B_function_engine; /// B-functions evaluator

    /// Calculate the Overlap integral between the two given STO basis functions
    /// \return the resulting energy quantity
    virtual energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &) ;

    /// Calculate the Overlap integral between the two given STO basis functions
    /// \return the resulting energy quantity
    virtual energy_unit_t kinetic(const std::array<STO_Basis_Function, 2> &) ;



    /// Given two specific B functions, calculate the overlap integral
    /// \param f1 First B function
    /// \param f2 Second B function
    /// \return Partial overlap integral value
    [[nodiscard]] std::complex<double> integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2)  override;

    /// This function is called back from the Gaussian Quadrature mechanism to get one value, at s, of the overlap
    /// integral.
    /// \param f1 details of the first B functions
    /// \param f2 details of the second B functions
    /// \param s point at which to calculate the integral
    /// \return overlap integral value at s
    static std::complex<double> calculate_overlap_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) ;



    static std::complex<double> get_gaunt_sum(const B_function_details &f1, const B_function_details &f2, double alpha) ;
    static std::complex<double> get_B_function_sum(const B_function_details &f1, const B_function_details &f2, double alpha, int l) ;
    static double calculate_delta(const B_function_details &f1, const B_function_details &f2, double s) ;
    static int get_l_min( const Quantum_Numbers &q1, const Quantum_Numbers &q2) ;
    static double get_gaunt_coeff(const std::array<const int, 6> &) ;
    static std::complex<double> calculate_B_function_value(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &point);

};

}