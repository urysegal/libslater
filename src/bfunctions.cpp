#include "bfunctions.h"
#include "logger.h"


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

double B_functions_representation_of_STO::calculate_coefficient(const Quantum_Numbers quantum_numbers,const unsigned int p) const
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
    auto Plm = bm::legendre_p(l,abs(m),std::cos(theta)); //Associated Legendre Polynomial


    auto Y = pow(i,m+abs(m));
    Y *= pow( ( (2*l + 1) * bm::factorial<double>(l-abs(m)) ) / (4*pi*(bm::factorial<double>(l+abs(m)))) , 1/2);
    Y *= Plm;
    Y *= std::exp(i* std::complex<double>(m*phi,1));

    return Y;
}//eval_spherical_harmonics

std::vector<double> B_function_Engine::cartesian_to_spherical(const center_t &r) const{
    bg::model::point<double, 3, bg::cs::cartesian> r_cart(r[0],r[1],r[2]);
    bg::model::point<double, 3, bg::cs::spherical<bg::radian>> r_spherical;
    bg::transform(r_cart, r_spherical);

    //RADIUS IS THIRD COORDINATE IN BOOST --(theta, phi, r)
    //We rearrange it as (r,theta,phi)
    std::vector<double> spher{r_spherical.get<2>(),r_spherical.get<0>(),r_spherical.get<1>()};

    return spher;
}//cartesian_to_spherical

std::complex<double> B_function_Engine::calculate(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &r) const
{
    //COMPLEX ARITHMETIC NEEDS TO BE Checked
    auto pi = bm::constants::pi<double>();
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;

    // Cartesian Representation of r to Spherical representation

    std::vector<double> spherical_coords = cartesian_to_spherical(r);
    double radius = spherical_coords[0];
    double phi = spherical_coords[1];
    double theta = spherical_coords[2];

    auto prefactor1 = pow(2.0/pi,1.0/2.0);
    auto prefactor2 = 1 / (pow(2.0,n+l) * bm::factorial<double>(n+l) );
    auto prefactor3 = pow(alpha*radius,(l+n-1.0/2.0)); //alpha*r needs to be corrected
    auto K = bm::cyl_bessel_k(n-1.0/2.0,alpha*radius); // May need to replace with recursion formula

    //need to extract theta and phi from r_spherical
    auto Y = eval_spherical_harmonics(quantum_numbers,theta,phi);

    return prefactor1*prefactor2*prefactor3 * K * Y;
}//B_function_Engine::calculate


}//namespace slater
