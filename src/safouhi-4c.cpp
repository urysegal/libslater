#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions.hpp>

#include "safouhi-4c.h"
#include "f2c.h"

using namespace std;
namespace bm = boost::math;

extern "C" {
int fortran_fourc(integer np1, integer l1, integer m1, doublereal zeta1,
                  integer np2, integer l2, integer m2, doublereal zeta2,
                  integer np3, integer l3, integer m3, doublereal zeta3,
                  integer np4, integer l4, integer m4, doublereal zeta4,
                  doublereal xa, doublereal ya, doublereal za,
                  doublereal xb, doublereal yb, doublereal zb,
                  doublereal xc, doublereal yc, doublereal zc,
                  doublereal xd, doublereal yd, doublereal zd,
                  doublereal *sdbar_res, doublereal *wgrep_res /* out vars */
);

double boost_choose(unsigned n, unsigned k)
{
    return bm::binomial_coefficient<double>(n,k)  ;
}

}

namespace slater
{

static Safouhi_4C_evaluator dummy_3c(safouhi_4c_name);

Safouhi_4C_evaluator::Safouhi_4C_evaluator() : STO_Integrator(4, 2, 0) {}

energy_unit_t
Safouhi_4C_evaluator::integrate(const std::vector<STO_Basis_Function> &functions, const std::vector<center_t> &centers)
{
    doublereal sdbar_res = 0;
    doublereal wgrep_res = 0;

    fortran_fourc(
        functions[0].get_quantum_numbers().n,
        functions[0].get_quantum_numbers().l,
        functions[0].get_quantum_numbers().m,
        functions[0].get_exponent(),

        functions[1].get_quantum_numbers().n,
        functions[1].get_quantum_numbers().l,
        functions[1].get_quantum_numbers().m,
        functions[1].get_exponent(),

        functions[2].get_quantum_numbers().n,
        functions[2].get_quantum_numbers().l,
        functions[2].get_quantum_numbers().m,
        functions[2].get_exponent(),

        functions[3].get_quantum_numbers().n,
        functions[3].get_quantum_numbers().l,
        functions[3].get_quantum_numbers().m,
        functions[3].get_exponent(),


        functions[0].get_center()[0],
        functions[0].get_center()[1],
        functions[0].get_center()[2],

        functions[1].get_center()[0],
        functions[1].get_center()[1],
        functions[1].get_center()[2],

        functions[2].get_center()[0],
        functions[2].get_center()[1],
        functions[2].get_center()[2],

        functions[3].get_center()[0],
        functions[3].get_center()[1],
        functions[3].get_center()[2],

        &sdbar_res, &wgrep_res /* out vars */
    );
    if ( this->use_sdbar ) {
        return sdbar_res;
    }
    return wgrep_res;
}


void Safouhi_4C_evaluator::init(const slater::STO_Integration_Options &params)
{
    std::string alg_name;
    params.get(Quadrature_4C_algorithm_Parameter_Name, alg_name);

    if ( alg_name == "WGREP") {
        this->use_sdbar = false;
    }

}

STO_Integrator *Safouhi_4C_evaluator::clone() const {
    return new Safouhi_4C_evaluator();
}

std::complex<double> Safouhi_4C_evaluator::integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2)
{
    // This is only relevant for 1-e integrals
    assert(false);
}

}