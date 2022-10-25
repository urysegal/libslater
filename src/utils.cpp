#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>

#include "slater-utils.h"
#include "libslater.h"

namespace bg = boost::geometry;
namespace bm = boost::math;

namespace slater
{

double distance(const center_t &A, const center_t &B)
{
    bg::model::point<double, 3, bg::cs::cartesian> A_point(A[0], A[1], A[2]);
    bg::model::point<double, 3, bg::cs::cartesian> B_point(B[0], B[1], B[2]);
    return bg::distance(A_point,B_point);
}

center_t vector_between(const center_t &A, const center_t &B)
{
    bg::model::point<double, 3, bg::cs::cartesian> A_point(A[0], A[1], A[2]);
    bg::model::point<double, 3, bg::cs::cartesian> B_point(B[0], B[1], B[2]);

    bg::subtract_point(B_point, A_point);

    return center_t({B_point.get<0>(), B_point.get<1>(), B_point.get<2>() });
}

center_t scale_vector( const center_t &P, double factor)
{
    return center_t {factor*P[0], factor*P[1], factor*P[2] } ;
}


std::complex<double> eval_spherical_harmonics(const Quantum_Numbers &quantumNumbers,const double theta,const double phi)
{
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
}
double compute_reduced_bessel_function_half(const double order,const double z){
    //auto pi = bm::constants::pi<double>();
    double k_tilde=0;
    auto n = int(order + 1.0/2.0);

    double k_tilde_0 = 0.5 * std::exp(-z);
    double k_tilde_1 = (1+z)/8 * std::exp(-z) ;

    // k_hat =(2/pi)^{1/2} *(alpha*r)^{n-1/2}* K_{n-1/2}(alpha*r)
    double k_hat=0;
    //if(z==0) {
        if (n == 1) {
            k_tilde = k_tilde_0;
        } else if (n == 2) {
            k_tilde = k_tilde_1;
        } else {
            for (int i = 3; i <= n; i++) {
                k_tilde = ((2.0 * double(i) - 3) / (2 * i)) * k_tilde_1 + ((z * z)/ (4.0 * double(i) * (i - 1)) * k_tilde_0);
                k_tilde_0 = k_tilde_1;
                k_tilde_1 = k_tilde;
            }
        }
         k_hat = pow(2, n) * bm::factorial<double>(n) * k_tilde;
    //} else {
      //   k_hat = pow(2.0/pi,1.0/2.0);
        // k_hat *= pow(z,( n - (1.0/2.0)));
         //k_hat *= bm::cyl_bessel_k(n-1.0/2.0,z);
    //}
    return k_hat;
}


double do_compute_reduced_bessel_function_half(
        const double order,
        const double z
)
{
    return compute_reduced_bessel_function_half(order, z);
}


}
