#include "bfunctions.h"
#include "logger.h"
#include "coordinates.h"
#include "slater-utils.h"

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
    double denominator =  bm::factorial<double>(2.0*p-n+l) * bm::double_factorial<double>(2.0*n-2.0*l-2.0*p);
    
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
    b_functions_sum_rescaling = std::pow(sto.get_exponent(), -int(quantum_numbers.n)+1);
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




std::complex<double> B_function_Engine::calculate(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &r) const
{
    // B_{n,l}^m(alpha,rr) = prefactor * K * Y
    // rr = (r,theta,phi)
    // prefactor =  1 / ( 2^{n+l} * (n+l)! )
    // k =(2/pi)^{1/2} *(alpha*r)^{n+1/2}* K_{n-1/2}(alpha*r)
    // Y = (alpha*r)^{l} * Y_{l,m}(theta,phi)

    //COMPLEX ARITHMETIC NEEDS TO BE Checked
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;

    // Cartesian Representation of r to Spherical representation
    Spherical_Coordinates spherical(r);

    // prefactor =  1 / ( 2^{n+l} * (n+l)! )
    auto prefactor = 1.0 / (pow(2.0,n+l) * bm::factorial<double>(n+l) );


     // May need to replace with recursion formula
    auto k_hat = compute_reduced_bessel_function_half(n-0.5,alpha*spherical.radius);
    // Y = (alpha*r)^{l} * Y_{l,m}(theta,phi)
    auto Y = eval_spherical_harmonics(quantum_numbers,spherical.theta,spherical.phi);
    Y *= pow(alpha*spherical.radius,(l));

    return prefactor * k_hat * Y;
}//calculate

Spherical_Coordinates::Spherical_Coordinates(const center_t &cartesian)
{
    auto x = cartesian[0];
    auto y = cartesian[1];
    auto z = cartesian[2];
    auto pi = bm::constants::pi<double>();

    radius = sqrt(x*x + y*y + z*z);
    if (radius == 0){
        theta=0;
        phi =0;
    }
    else {
        theta = acos(z / radius);
        if (x > 0) {
            phi = atan(y / x);
        } else if (x < 0) {
            if (y >= 0) {
                phi = atan(y / x) + pi;
            } else {
                phi = atan(y / x) - pi;
            }
        }
            // x==0
        else {
            if (y > 0) {
                phi = pi / 2;
            } else if (y < 0) {
                phi = -pi / 2;
            }
                //x==0, y==0, phi = undefined
            else {
                phi = 0;
            }
        }
    }
}


}//namespace slater
