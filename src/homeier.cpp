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
        energy_unit_t partial_result = integrate_overlap_using_b_functions(p.first, p.second);
        partial_results.emplace_back(partial_result);
    }

    energy_unit_t final_result = 0 ;
    for ( auto &pr : partial_results )
        final_result += pr;

    return final_result;
}

energy_unit_t Homeier_Integrator::integrate_overlap_using_b_functions(const B_function_details &f1, const B_function_details &f2) const
{
    auto f = [&](const double& s) { return this->calculate_overlap_gaussian_point(f1, f2, s) ;};
    double Q = boost::math::quadrature::gauss<double, 7>::integrate(f, 0, 1);
    return Q;
}


double Homeier_Integrator::calculate_overlap_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const
{
    double W_hat = calculate_W_hat(f1, f2, s);
    double S = calculate_S(f1, f2 ,s);
    return W_hat * S;
}



double Homeier_Integrator::calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s) const
{
    //TESTING NEEDED
    //Equation 30 in On the Evaluation of Overlap Integrals with Exponential-Type Basis Functions
    //HERBERT H. H. HOMEIER AND E. OTTO STEINBORN
    //W_hat(s) is the new weight function after mobius transformation
    auto quantum_numbers_1 = f1.get_quantum_numbers();
    auto quantum_numbers_2 = f2.get_quantum_numbers();

    auto n1 = quantum_numbers_1.n;
    auto l1 = quantum_numbers_1.l;
    
    auto n2 = quantum_numbers_2.n;
    auto l2 = quantum_numbers_2.l;

    auto alpha = f1.get_exponent();
    auto beta = f2.get_exponent();
    
    double eta = beta/alpha;

    double numerator = pow(1-s,n1+l1) * pow(s,n2+l2);
    double denominator1 = pow(s+(1-s)*eta,(l1+l2+3)/2);
    double denominator2 = pow( (1-s)*eta*alpha*alpha + s*beta*beta,n1+n2+(l1+l2+1)/2);
    double prefactor  = pow(eta, n1+l1+1);
    
    return prefactor * numerator / (denominator1*denominator2);
}

double Homeier_Integrator::calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const
{
    return f2.get_quantum_numbers().n * f1.get_coefficient() + s;
}



}