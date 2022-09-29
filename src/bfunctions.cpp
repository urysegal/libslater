#include "bfunctions.h"
#include "logger.h"
#include <iostream> //REMOVE - only added for my primitive caveman style debugging
#include <boost/math/special_functions/factorials.hpp>

namespace bmath = boost::math;

namespace slater
{


B_function_details::B_function_details(sto_coefficient_t coefficient_, sto_exponent_t exponent_,
                                       const Quantum_Numbers &quantum_numbers_, const center_t &center_) :
                                       coefficient(coefficient_), exponent(exponent_),
                                       quantum_numbers(quantum_numbers_), center(center_)
{}//B_function_details


sto_coefficient_t B_function_details::get_coefficient() const
{
    return coefficient;
}

sto_exponent_t B_function_details::get_exponent() const
{
    return exponent;
}

const Quantum_Numbers &B_function_details::get_quantum_numbers() const
{
    return quantum_numbers;
}

double B_functions_representation_of_STO::calculate_coefficient(const B_function_details &bfd,const unsigned int p) const
{
    //TESTING NEEDED

    // B function representation of an STO is 
    // X_{n,l}^m (alpha,r) = N * sum_p ^{n-l} Bcoeff * B_{p,l}^m (alpha,r) // MAY CHANGE REPRESENTATION 
    //calculate the individual coefficients Bcoeff

     
    auto quantum_numbers = bfd.get_quantum_numbers();
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;

    //the factorial function call is too verbose - shorten it ?
    double numerator = pow(-1.0,(n-l-p)) * bmath::factorial<double>(n-l) * pow(2.0,(l+p)) * bmath::factorial<double>(l+p) ;
    double denominator =  bmath::factorial<double>(2.0*p-n+l) * bmath::factorial<double>((boost::math::factorial<double>(2.0*n-2.0*l-2.0*p)));
    
    return numerator/denominator;
}//calculate_coefficient

B_functions_representation_of_STO::B_functions_representation_of_STO(const STO_Basis_Function &sto)
{
    //TESTING NEEDED
    // B function representation of an STO is 
    // X_{n,l}^m (alpha,r) = N * sum_p ^{n-l} Bcoeff * B_{p,l}^m (alpha,r)
    // p depends on even/odd n-l
    // Bcoeff is computed via calculate_coefficient  
    auto quantum_numbers = sto.get_quantum_numbers();
    auto n = quantum_numbers.n;
    auto l = quantum_numbers.l;
    unsigned int min_p;

    if ( (n-l)%2 ==0 ){
        min_p = (n-l)/2; 
    }
    else{
        min_p = (n-l+1)/2;
    }

    //loop over p to create all Bcoeffs and Bfunctions
    for (unsigned int p= min_p; p<= (n-l) ; p++){
     B_function_details B_func_p = B_function_details(sto.get_coefficient(), sto.get_exponent()/2, quantum_numbers, sto.get_center());
     double Bcoeff_p = calculate_coefficient(B_func_p,p);

     //store in vectors 
     components.emplace_back(B_func_p);
     components_coefficients.emplace_back(Bcoeff_p);
    }

}//B_functions_representation_of_STO

}//namespace slater
