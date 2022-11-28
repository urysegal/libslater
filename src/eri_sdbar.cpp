#include "eri_sdbar.h"


namespace slater {

static Electron_Repulsion_SDbar dummy_eri_sdbar(electron_repulsion_sdbar_name);


Electron_Repulsion_SDbar::Electron_Repulsion_SDbar() : STO_Integrator(3, 1, 1) {}

void Electron_Repulsion_SDbar::init(const slater::STO_Integration_Options &params) {
    params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
}

STO_Integrator *Electron_Repulsion_SDbar::clone() const {
    return new Electron_Repulsion_SDbar();
}


energy_unit_t Electron_Repulsion_SDbar::integrate(const std::vector<STO_Basis_Function> &functions,
                                                  const std::vector<center_t> &centers)
{
    setup_state(functions);
    return evaluate();
}


void Electron_Repulsion_SDbar::setup_state(const std::vector<STO_Basis_Function> &functions)
{

    state.A = functions[0].get_center();
    state.B = functions[1].get_center();
    state.C = functions[2].get_center();
    state.D = functions[3].get_center();
    for ( auto i = 0U ; i < 4 ; ++i ) {
        state.zeta[i] = functions[i].get_exponent();
        state.n[i] = functions[i].get_quantum_numbers().n;
        state.l[i] = functions[i].get_quantum_numbers().l;
        state.m[i] = functions[i].get_quantum_numbers().m;
    }

}


}