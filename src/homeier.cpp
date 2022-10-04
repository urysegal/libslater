#include <boost/math/quadrature/gauss.hpp>


#include "homeier.h"
#include "bfunctions.h"
#include "logger.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/binomial.hpp>

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
        energy_unit_t partial_result = integrate_overlap_using_b_functions(p.first.second, p.second.second);
        partial_result *= p.first.first*p.second.first;
        partial_results.emplace_back(partial_result);
    }

    energy_unit_t final_result = 0 ;
    for ( auto &pr : partial_results )
        final_result += pr;

    final_result *= functions[0].get_coefficient() * functions[1].get_coefficient() ;

    final_result *= f1.get_rescaling_coefficient() * f2.get_rescaling_coefficient() ;

    return final_result;
}

std::complex<double> Homeier_Integrator::integrate_overlap_using_b_functions(const B_function_details &f1, const B_function_details &f2) const
{
    auto f = [&](const double& s) { return this->calculate_overlap_gaussian_point(f1, f2, s) ;};
    std::complex<double> Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);
    //need to add multiply with another prefactor here - Equation 20
    return Q;
}


    std::complex<double> Homeier_Integrator::calculate_overlap_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const
{
    double W_hat = calculate_W_hat(f1, f2, s);
    std::complex<double> S = calculate_S(f1, f2 ,s);
    return std::complex<double>(W_hat,0) * S;
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

    auto alpha = f1.get_alpha();
    auto beta = f2.get_alpha();
    
    double eta = beta/alpha;

    double numerator = pow(1-s,n1+l1) * pow(s,n2+l2);
    double denominator1 = pow(s+(1-s)*eta,(l1+l2+3)/2);
    double denominator2 = pow( (1-s)*eta*alpha*alpha + s*beta*beta,n1+n2+(l1+l2+1)/2);
    double prefactor  = pow(eta, n1+l1+1);
    
    return prefactor * numerator / (denominator1*denominator2);
}

double Homeier_Integrator::calculate_delta(const B_function_details &f1, const B_function_details &f2, double s) const
{
    return f1.get_alpha() * f2.get_quantum_numbers().l * s;
}

std::complex<double> Homeier_Integrator::get_B_function_sum(const B_function_details &f1, const B_function_details &f2, double alpha, int l) const
{
    std::complex<double> total_sum = 0;
    const Quantum_Numbers &q1 = f1.get_quantum_numbers();
    const Quantum_Numbers &q2 = f2.get_quantum_numbers();

    int delta_l = (q1.l + q2.l - l );
    assert( delta_l % 2 == 0);
    delta_l/=2;

    for (auto j = 0 ; j < delta_l ; j++) {
        Quantum_Numbers B_function_parameters = {q1.n+q2.n+2*delta_l+1-j, (unsigned int)l, q2.m - q1.m };
        total_sum += std::complex<double>(pow(-1, j) * boost::math::binomial_coefficient<double>(delta_l, j),0) *
                calculate_B_function_value(B_function_parameters, alpha, f2.get_center());
    }
    return total_sum; }

std::complex<double> Homeier_Integrator::get_gaunt_sum(const B_function_details &f1, const B_function_details &f2, double alpha) const
{

    const Quantum_Numbers &q1 = f1.get_quantum_numbers();
    const Quantum_Numbers &q2 = f2.get_quantum_numbers();

    int l_min = get_l_min(q1, q2);

    int l_max = f1.get_quantum_numbers().l + f2.get_quantum_numbers().l;

    std::complex<double> total_sum = 0;
    for ( auto l = l_min ; l < l_max ; l+=2 )
    {
        auto gaunt_coeff = get_gaunt_coeff({(int)q2.l, q2.m, (int)q1.l,q1.m, l, q2.m - q1.m});
        std::complex<double> B_function_sum = get_B_function_sum(f1, f2, alpha, l);
        total_sum += std::complex<double>(gaunt_coeff,0) * B_function_sum;
    }
    return total_sum;
}


std::complex<double> Homeier_Integrator::calculate_S(const B_function_details &f1, const B_function_details &f2, double s) const
{
    // formula 13 second paper or 19 first paper
    // We compute S_{n1,l1,m1}^{n2,l2,m2}(d,d,R) using Equation 13
    // d = delta(alpha,beta,s)
    // We will need Gaunt coefficients here

    auto l2 = f2.get_quantum_numbers().l;
    auto delta = calculate_delta(f1,f2,s);
    auto pi = boost::math::constants::pi<double>();


    auto gaunt_sum = get_gaunt_sum(f1, f2, delta);


    std::complex<double> result =std::complex<double>(pow(-1, l2 ) * (4*pi/pow(delta,3)),0) * gaunt_sum ;

    return result;
}

double Homeier_Integrator::get_gaunt_coeff(const std::array<const int, 6> &args) const
{
    return gaunt_engine.calculate(args);
}


std::complex<double> Homeier_Integrator::calculate_B_function_value(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &point) const
{
    return B_function_engine.calculate(quantum_numbers, alpha, point);
}


int Homeier_Integrator::get_l_min( const Quantum_Numbers &q1, const Quantum_Numbers &q2) const
{
    return std::min(q1.l, q2.l); // Wrong , Gautam please write
}






}