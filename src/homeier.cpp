#include <boost/math/quadrature/gauss.hpp>


#include "homeier.h"
#include "bfunctions.h"
#include "logger.h"

namespace slater {

Homeier_Integrator::Homeier_Integrator() : STO_Integrator()
{

}

Homeier_Integrator::~Homeier_Integrator()
{
}

void Homeier_Integrator::init(const STO_Integration_Options &params)
{
    params.get(Use_Normalized_B_Functions_Parameter_Name, use_normalized_b_functions);
    params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
}


void Homeier_Integrator::create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2)
{
    for ( auto i : f1) {
        for (auto j: f2) {
            equivalence_series.emplace_back(i, j);
        }
    }
}

energy_unit_t Homeier_Integrator::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    B_functions_representation_of_STO f1(functions[0]);
    B_functions_representation_of_STO f2(functions[1]);

    create_integration_pairs(f1, f2);

    std::vector<energy_unit_t> partial_results;
    for ( auto const &p : equivalence_series )
    {
        energy_unit_t partial_result = integrate_with_b_functions(p.first, p.second);
        partial_results.emplace_back(partial_result);
    }

    energy_unit_t final_result = 0 ;
    for ( auto &pr : partial_results )
        final_result += pr;

    return final_result;
}

energy_unit_t Homeier_Integrator::integrate_with_b_functions(const B_function_details &f1, const B_function_details &f2) const
{
    auto f = [&](const double& s) { return this->calculate_gaussian_point(f1, f2, s) ;};
    double Q = boost::math::quadrature::gauss<double, 7>::integrate(f, 0, 1);
    return Q;
}


double Homeier_Integrator::calculate_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const
{
    double W_hat = calculate_W_hat(f1, f2, s);
    double S = calculate_S(f1, f2 ,s);
    return W_hat * S;
}



double Homeier_Integrator::calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s) const
{
    /// Gautam - Equation 30 in the second paper
    return f1.get_exponent() * f2.get_quantum_numbers().l * s;
}

double Homeier_Integrator::calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const
{
    return f2.get_quantum_numbers().n * f1.get_coefficient() + s;
}



}