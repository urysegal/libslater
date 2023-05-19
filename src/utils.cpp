#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>
#include "coordinates.h"
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

double vector_length(const center_t &A)
{
    auto mag = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    return mag;
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

double pochhammer(const double x, const double n){
    return std::tgamma(x+n) / std::tgamma(x);
}
std::complex<double> eval_spherical_harmonics(const Quantum_Numbers &quantumNumbers,const Spherical_Coordinates &spherical)
{
    // Evaluates Spherical Harmonics Y_l^m(theta,phi)
    auto m = quantumNumbers.m;
    auto l = quantumNumbers.l;

//   compute spherical harmonic for m
    auto Y = bm::spherical_harmonic(l,m,spherical.theta,spherical.phi);

    return Y;
}

    std::complex<double> eval_spherical_harmonics_real(const Quantum_Numbers &quantumNumbers,const Spherical_Coordinates &spherical)
    {
        // Evaluates Spherical Harmonics Y_l^m(theta,phi)
        auto m = quantumNumbers.m;
        auto l = quantumNumbers.l;

//   compute spherical harmonic for POSITIVE m
        auto Y = bm::spherical_harmonic(l,m,spherical.theta,spherical.phi);

//   compute spherical harmonic for NEGATIVE m
        auto Yn = bm::spherical_harmonic(l,-1*m,spherical.theta,spherical.phi);

        m = quantumNumbers.m;
        if (m<0){
            // Ylm  = i/sqrt(2) * (Yl(m) - (-1)^m * Yl(-m))
            Y = std::complex<double>(0,1)*(1.0/sqrt(2.0))*(Y-pow(-1,m) *Yn);
        }
        else if (m>0){
            // Ylm  = 1/sqrt(2) * (Yl(-m) + (-1)^m * Yl(m))
            Y = (1.0/sqrt(2.0))*(Yn + pow(-1,m) *Y);
        }
        return Y;
    }

double compute_reduced_bessel_function_half(const double order,const double z){
    auto pi = bm::constants::pi<double>();
    double k_tilde=0;
    auto n = int(order + 1.0/2.0);

    double k_tilde_0 = 0.5 * std::exp(-z);
    double k_tilde_1 = (1+z)/8 * std::exp(-z) ;

    //k_hat =(2/pi)^{1/2} *(alpha*r)^{n-1/2}* K_{n-1/2}(alpha*r)
    double k_hat=0;
    if(z==0) {
        if (n == 1) {
            k_tilde = k_tilde_0;
        } else if (n == 2) {
            k_tilde = k_tilde_1;
        } else {
            for (int i = 3; i <= n; i++) {
                k_tilde = ((2.0 * double(i) - 3.0) / (2.0 * i)) * k_tilde_1 + ((z * z)/ (4.0 * double(i) * (i - 1.0)) * k_tilde_0);
                k_tilde_0 = k_tilde_1;
                k_tilde_1 = k_tilde;
            }
        }
         k_hat = pow(2, n) * bm::factorial<double>(n) * k_tilde;
    } else {
         k_hat = pow(2.0/pi,1.0/2.0);
         k_hat *= pow(z,( n - (1.0/2.0)));
         k_hat *= bm::cyl_bessel_k(n-1.0/2.0,z);
    }
    return k_hat;
}


}
