#include <boost/math/quadrature/gauss.hpp>


#include "homeier.h"
#include "bfunctions.h"
#include "logger.h"

namespace slater {

Homeier_Integrator::~Homeier_Integrator()
{

}

void Homeier_Integrator::init(const STO_Integration_Options &params)
{
    params.get(Use_Normalized_B_Functions_Parameter_Name, use_normalized_b_functions);
    params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
}


void Homeier_Integrator::create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2) const
{
    for ( int i = 0 ; i < f1.size() ; i++)
        for ( int j = 0 ; i < fj.size() ; i++ )
            components.emplace_back(f1.get(i), f2.get(j));

}

integral_value Homeier_Integrator::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    B_functions_representation_of_STO f1(functions[0]);
    B_functions_representation_of_STO f2(functions[1]);

    create_integration_pairs();

    std::vector<double> partial_results;
    for ( auto const &p : components )
    {
        partial_results.emplace_back(integrate_with_b_functions(p.first, p.second));
    }

    double final_result = 0 ;
    for ( auto &pr : partial_results )
        final_result += pr;

    return final_result;
}

double Homeier_Integrator::calculate_one_weight(const B_function_details &f1, const B_function_details &f2, double s)
{
    double W_hat = calculate_W_hat(f1, f2, s);
    double S = calculate_S(f1, f2 ,s);
    return W_hat * S;
}


double Homeier_Integrator::integrate_with_b_functions(const B_function_details &f1, const B_function_details &f2)
{
    auto f = [&](const double& s) { this->calculate_one_weight(f1, f2, s) };
    double Q = boost::math::quadrature::gauss<double, 7>::integrate(f, 0, 1);
    return Q;
}


}