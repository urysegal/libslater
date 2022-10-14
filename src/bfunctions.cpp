#include "bfunctions.h"
#include "logger.h"
#include "coordinates.h"

namespace bm = boost::math;
namespace bg = boost::geometry;
namespace slater
{


B_function_details::B_function_details( sto_exponent_t exponent_,
                                       const Quantum_Numbers &quantum_numbers_, const center_t &center_) :
        alpha(exponent_),
        quantum_numbers(quantum_numbers_), center(center_)
{}//B_function_details



sto_exponent_t B_function_details::get_alpha() const
{
    return alpha;
}

const Quantum_Numbers &B_function_details::get_quantum_numbers() const
{
    return quantum_numbers;
}

double B_functions_representation_of_STO::calculate_coefficient(const Quantum_Numbers &quantum_numbers,const unsigned int p) const
{
    //TESTING NEEDED

    // B function representation of an STO is 
    // X_{n,l}^m (alpha,r) = N * sum_p ^{n-l} Bcoeff * B_{p,l}^m (alpha,r)
    //calculate the individual coefficients Bcoeff

    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;

    //the factorial function call is too verbose - shorten it ?
    double numerator = pow(-1.0,(n-l-p)) * bm::factorial<double>(n-l) * pow(2.0,(l+p)) * bm::factorial<double>(l+p) ;
    double denominator =  bm::factorial<double>(2.0*p-n+l) * bm::factorial<double>(bm::factorial<double>(2.0*n-2.0*l-2.0*p));
    
    return numerator/denominator;
}//calculate_coefficient

B_functions_representation_of_STO::B_functions_representation_of_STO(const STO_Basis_Function &sto, const center_t &new_center)
{
    //TESTING NEEDED
    // B function representation of an STO is 
    // X_{n,l}^m (alpha,r) = N * sum_p ^{n-l} Bcoeff * B_{p,l}^m (alpha,r)
    // p depends on even/odd n-l
    // Bcoeff is computed via calculate_coefficient
    auto quantum_numbers = sto.get_quantum_numbers();
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;
    b_functions_sum_rescaling = std::pow(sto.get_exponent(), (-quantum_numbers.n)+1);

    unsigned int min_p;

    if ( (n-l)%2 ==0 ){
        min_p = (n-l)/2; 
    }
    else{
        min_p = (n-l+1)/2;
    }

    //loop over p to create all Bcoeffs and Bfunctions
    auto quantum_numbers_p = sto.get_quantum_numbers();
    for (unsigned int p = min_p; p <= (n - l); p++) {
        //Create  B_{p,l}^m (alpha,r)
        quantum_numbers_p.n = p;
        B_function_details B_func_p = B_function_details(sto.get_exponent(), quantum_numbers_p, new_center);

        //Compute summing coefficient for B_{p,l}^m (alpha,r)
        double Bcoeff_p = calculate_coefficient(quantum_numbers, p);

        //store in vectors
        components.emplace_back(Bcoeff_p,B_func_p);
    }

}//B_functions_representation_of_STO

std::complex<double> B_function_Engine::eval_spherical_harmonics(const Quantum_Numbers quantumNumbers,const double theta,const double phi) const{
    // Evaluates Spherical Harmonics Y_l^m(theta,phi)
    auto pi = bm::constants::pi<double>();
    std::complex<double>  i(0,1);

    auto m = quantumNumbers.m;
    auto l = quantumNumbers.l;
    auto Plm = bm::legendre_p(l,m,std::cos(theta)); //Associated Legendre Polynomial

    //Using Wikipedia's accoustics definition to stay consistent with legendre_p function in boost
    // which includes the Condon-Shortley phase term
    std::complex<double> Y;
    Y = pow( ( (2*l + 1) * bm::factorial<double>(l-m) ) / (4 * pi * bm::factorial<double>(l + m)) , 1.0 / 2);
    Y *= Plm;
    Y *= std::exp(std::complex<double>(0,m*phi));

    return Y;
}//eval_spherical_harmonics



std::complex<double> B_function_Engine::calculate(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &r) const
{
    // B_{n,l}^m(alpha,rr) = prefactor * K * Y
    // rr = (r,theta,phi)
    // prefactor = (2/pi)^{1/2} * (2^{n+l} * (n+l)! )^{-1} (alpha*r)^{l+n+1/2}
    // K = K_{n-1/2}(alpha*r)
    // Y = Y_{l,m}(theta,phi)

    //COMPLEX ARITHMETIC NEEDS TO BE Checked
    auto pi = bm::constants::pi<double>();
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;

    // Cartesian Representation of r to Spherical representation
    Spherical_Coordinates spherical(r);

    auto prefactor = pow(2.0/pi,1.0/2.0);
    prefactor *= 1 / (pow(2.0,n+l) * bm::factorial<double>(n+l) );
    prefactor *= pow(alpha*spherical.radius,(l+n-1.0/2.0)); //alpha*r needs to be corrected

    //modified Bessel Function of Second Kind
    auto K = bm::cyl_bessel_k(n-1.0/2.0,alpha*spherical.radius); // May need to replace with recursion formula

    //need to extract theta and phi from r_spherical
    auto Y = eval_spherical_harmonics(quantum_numbers,spherical.theta,spherical.phi);

    return prefactor * K * Y;
}//calculate

Spherical_Coordinates::Spherical_Coordinates(const center_t &cartesian)
{
    bg::model::point<double, 3, bg::cs::cartesian> r_cart(cartesian[0],cartesian[1],cartesian[2]);
    bg::model::point<double, 3, bg::cs::spherical<bg::radian>> r_spherical;
    r_spherical.set<1>(0);
    bg::transform(r_cart, r_spherical);

    theta = r_spherical.get<0>();
    phi = r_spherical.get<1>();
    radius = r_spherical.get<2>();
}


}//namespace slater
