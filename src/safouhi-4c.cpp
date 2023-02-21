#include "safouhi-4c.h"
#include "f2c.h"

using namespace std;

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
}

namespace slater
{

static Safouhi_4C_evaluator dummy_3c(safouhi_4c_name);

Safouhi_4C_evaluator::Safouhi_4C_evaluator() : STO_Integrator(4, 2, 0) {}

void Safouhi_4C_evaluator::init(const slater::STO_Integration_Options &params)
{
    std::string alg_name;
    params.get(Quadrature_4C_algorithm_Parameter_Name, alg_name);

    if ( alg_name == "GREPW") {
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