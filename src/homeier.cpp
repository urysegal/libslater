#include <boost/math/quadrature/gauss.hpp>


#include "homeier.h"
#include "bfunctions.h"
#include "logger.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace bm = boost::math;
namespace slater {

Homeier_Integrator dummy_integrator(overlap_homeier_imp_name);

B_function_Engine Homeier_Integrator::B_function_engine;

STO_Integrator *Homeier_Integrator::clone() const
{
    return new Homeier_Integrator();
}


Homeier_Integrator::~Homeier_Integrator()
{
}

energy_unit_t Homeier_Integrator::integrate( const std::vector<STO_Basis_Function> &functions, const std::vector<center_t> &centers)
{
    return overlap({functions[0], functions[1]});
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
    /// R = R_2 - R_1
    new_centers[1] = {c2[0]-c1[0],c2[1]-c1[1],c2[2]-c1[2]};
    new_centers[0] = {0,0,0};
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
        // int(B_1*B_2)
        energy_unit_t partial_result = integrate_overlap_using_b_functions(p.first.second, p.second.second);

        // Bcoeff1*Bcoeff2 * int(B_1*B_2)
        partial_result *= p.first.first*p.second.first;
        partial_results.emplace_back(partial_result);
    }

    // sum(sum(Bcoeff1*Bcoeff2 * int(B_1 B_2)))
    energy_unit_t final_result = 0 ;
    for ( auto &pr : partial_results )
        final_result += pr;

    // alpha_1^{-n1+1} * alpha_2^{-n2+1} sum(sum(Bcoeff1*Bcoeff2 * int(B_1 B_2)))
    final_result *= functions[0].get_coefficient() * functions[1].get_coefficient() ;

    // N1*N2 * alpha_1^{-n+1} * alpha_2^{-n+1} sum(sum(Bcoeff1*Bcoeff2 * int(B_1 B_2)))
    final_result *= f1.get_rescaling_coefficient() * f2.get_rescaling_coefficient() ;

    return final_result;
}

std::complex<double> Homeier_Integrator::integrate_overlap_using_b_functions(const B_function_details &f1, const B_function_details &f2) const
{
    auto f = [&](const double& s) { return this->calculate_overlap_gaussian_point(f1, f2, s) ;};
    // int(W_hat*S)
    std::complex<double> Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);

    // prefactor*int(W_hat*S)
    auto q1 = f1.get_quantum_numbers();
    auto q2 = f2.get_quantum_numbers();
    auto alpha = f1.get_alpha();
    auto beta  = f2.get_alpha();
    double prefactor = pow(alpha,2*q1.n+q1.l-1)*pow(beta,2*q2.n+q2.l-1);
    prefactor *= bm::factorial<double>(q1.n + q1.l + q2.n + q2.l + 1 );
    prefactor /=  bm::factorial<double>(q1.n + q1.l) * bm::factorial<double>(q2.n + q2.l);
    Q *= std::complex<double>(prefactor,0);
    return Q;
}


    std::complex<double> Homeier_Integrator::calculate_overlap_gaussian_point(const B_function_details &f1, const B_function_details &f2, double s) const
{
    double W_hat = calculate_W_hat(f1, f2, s);
    std::complex<double> S = calculate_S(f1, f2 ,s);
    return std::complex<double>(W_hat,0) * S;
}



double Homeier_Integrator::calculate_W_hat(const B_function_details &f1, const B_function_details &f2, double s)
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

    double numerator = pow(1.0-s,n1+l1) * pow(s,n2+l2);
    double denom_1_power = (3.0 + l1 + l2 ) / 2.0 ;
    double denominator1 = pow(s+(1.0-s)*eta, denom_1_power);
    double denom_2_power = n1+n2+double(l1+l2+1.0)/2.0 ;
    double denominator2 = pow( (1.0-s)*eta*alpha*alpha + s*beta*beta, denom_2_power);
    double prefactor  = pow(eta, n1+l1+1);
    
    return prefactor * numerator / (denominator1*denominator2);
}

double Homeier_Integrator::calculate_delta(const B_function_details &f1, const B_function_details &f2, double s)
{
    //delta(alpha,beta,t)
    auto alpha = f1.get_alpha();
    auto beta = f2.get_alpha();

    auto eta = beta / alpha ;
    double t = s/( s+ (1.0-s) * eta );

    auto delta = sqrt((1.0-t)*alpha*alpha + t*beta*beta);
    return delta;
}

std::complex<double> Homeier_Integrator::get_B_function_sum(const B_function_details &f1, const B_function_details &f2, double alpha, int l)
{
    std::complex<double> total_sum = 0;
    const Quantum_Numbers &q1 = f1.get_quantum_numbers();
    const Quantum_Numbers &q2 = f2.get_quantum_numbers();

    int delta_l = (q1.l + q2.l - l );

    assert( delta_l % 2 == 0);
    delta_l/=2;

    for (auto j = 0 ; j <= delta_l ; j++) {
        Quantum_Numbers B_function_parameters = {q1.n+q2.n+2*delta_l+1-j, l, q2.m - q1.m };
        total_sum += std::complex<double>(pow(-1, j) * boost::math::binomial_coefficient<double>(delta_l, j),0) *
                calculate_B_function_value(B_function_parameters, alpha, f2.get_center());
    }
    return total_sum;
}

std::complex<double> Homeier_Integrator::get_gaunt_sum(const B_function_details &f1, const B_function_details &f2, double alpha)
{

    const Quantum_Numbers &q1 = f1.get_quantum_numbers();
    const Quantum_Numbers &q2 = f2.get_quantum_numbers();

    int l_min = get_l_min(q1, q2);

    int l_max = f1.get_quantum_numbers().l + f2.get_quantum_numbers().l;

    std::complex<double> total_sum = 0;
    for ( auto l = l_min ; l <= l_max ; l+=2 )
    {
        auto gaunt_coeff = get_gaunt_coeff({(int)q2.l, q2.m, (int)q1.l,q1.m, l, q2.m - q1.m});
        std::complex<double> B_function_sum = get_B_function_sum(f1, f2, alpha, l);
        total_sum += gaunt_coeff * B_function_sum;
    }
    return total_sum;
}


std::complex<double> Homeier_Integrator::calculate_S(const B_function_details &f1, const B_function_details &f2, double s)
{
    // formula 13 second paper or 19 first paper
    // We compute S_{n1,l1,m1}^{n2,l2,m2}(d,d,R) using Equation 13
    // d = delta(alpha,beta,s)
    // We will need Gaunt coefficients here

    auto l2 = f2.get_quantum_numbers().l;
    auto delta = calculate_delta(f1,f2,s);
    auto pi = boost::math::constants::pi<double>();


    auto gaunt_sum = get_gaunt_sum(f1, f2, delta);


    std::complex<double> result = pow(-1, l2 ) * (4.0*pi) * gaunt_sum ;

    return result;
}

double Homeier_Integrator::get_gaunt_coeff(const std::array<const int, 6> &args)
{
    return Gaunt_Coefficient_Engine::get()->calculate(args);
}


std::complex<double> Homeier_Integrator::calculate_B_function_value(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &point)
{
    return B_function_engine.calculate(quantum_numbers, alpha, point);
}


int Homeier_Integrator::get_l_min( const Quantum_Numbers &q1, const Quantum_Numbers &q2)
{
    /// l_min depends on whether max( |l1 - l2|, |m1-m2|) + l_max is even or odd

    auto m = std::max(abs(int(q1.l-q2.l)),abs(q1.m-q2.m));
    auto switch_condition = m +q1.l+q2.l;
    return m + (switch_condition%2) ;
}


/// These functions are only for debug purpose


void calculate_gauss_point(const STO_Basis_Function &bf1, const STO_Basis_Function &bf2, double s)
{

    center_t new_centers[2];
    shift_first_center_to_origin(bf1.get_center(), bf2.get_center(), new_centers);

    B_functions_representation_of_STO f1(bf1, new_centers[0]);
    B_functions_representation_of_STO f2(bf2, new_centers[1]);


    for (auto i: f1) {
        for (auto j: f2) {
            std::complex<double> S = Homeier_Integrator::calculate_S(i.second, j.second, s);
            auto w_hat = Homeier_Integrator::calculate_W_hat(i.second, j.second, s);
            fprintf(stdout, "S= %15.15f + %15.15f i ,  W=%15.15f  \n", S.real(), S.imag(), w_hat );
        }
    }
}

}

