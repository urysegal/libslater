#include "bfunctions.h"
#include "logger.h"
#include <iostream> //REMOVE - only added for my primitive caveman style debugging
#include <boost/math/special_functions/factorials.hpp>

namespace bmath = boost::math;

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
    double numerator = pow(-1.0,(n-l-p)) * bmath::factorial<double>(n-l) * pow(2.0,(l+p)) * bmath::factorial<double>(l+p) ;
    double denominator =  bmath::factorial<double>(2.0*p-n+l) * bmath::factorial<double>((boost::math::factorial<double>(2.0*n-2.0*l-2.0*p)));
    
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


double B_function_Engine::calculate(const Quantum_Numbers &quantum_numbers, double alpha, const center_t &r) const
{
    return quantum_numbers.n + alpha + r[1]; /// Gotham please implement
}


}//namespace slater
