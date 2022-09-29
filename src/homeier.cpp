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

/// Shift the centers of the B functions so that one for f1 is (0,0,0) and the other shifted accordingly
/// \param c1  First center to shift. This one shifts to (0,0,0)
/// \param c2 Second center to shift
/// \param new_centers
void shift_first_center_to_origin(const center_t &c1, const center_t c2, center_t *new_centers)
{
    /// Gautam - make the right calculation
    new_centers[0] = c1;
    new_centers[1] = c2;
}


energy_unit_t Homeier_Integrator::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    center_t new_centers[2];
    shift_first_center_to_origin(functions[0].get_center(), functions[1].get_center(), new_centers);

    B_functions_representation_of_STO f1(functions[0] , new_centers[0]);
    B_functions_representation_of_STO f2(functions[1], new_centers[1]);

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

double Homeier_Integrator::integrate_overlap_using_b_functions(const B_function_details &f1, const B_function_details &f2) const
{
    auto f = [&](const double& s) { return this->calculate_overlap_gaussian_point(f1, f2, s) ;};
    double Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);
    //need to add multiply with another prefactor here - Equation 20
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
    //eta is computed as per Equation 35

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
    // We compute S_{n1,l1,m1}^{n2,l2,m2}(d,d,R) using Equation 13
    // d = delta(alpha,beta,s)
    // We will need Gaunt coefficients here

    return f2.get_quantum_numbers().n  + s;
}



}