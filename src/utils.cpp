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


}
