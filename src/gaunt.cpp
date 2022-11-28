#include <math.h>
#include "gaunt.h"
#include "../precalculate/gaunt-table.h"
#include <assert.h>

namespace slater {

Gaunt_Coefficient_Engine *Gaunt_Coefficient_Engine::instance = nullptr;


Gaunt_Coefficient_Engine *
Gaunt_Coefficient_Engine::get()
{
    if ( not instance ) {
        instance = new Gaunt_Coefficient_Engine();
    }
    return instance;
}


Gaunt_Coefficient_Engine::Gaunt_Coefficient_Engine()
{}

int Gaunt_Coefficient_Engine::get_maximal_gaunt_l()
{
    return  L_MAX;
}


double Gaunt_Coefficient_Engine::calculate(const std::array<const int, 6> &args) const
{
    auto l1 = args[0];

    auto orig_m1 = args[1];

    auto m1 = -orig_m1 + L_MAX;

    auto l2 = args[2];
    auto m2 = args[3] + L_MAX;

    auto l3 = args[4];
    auto m3 = args[5] + L_MAX;

    assert(l1>=0);
    assert(l2>=0);
    assert(l3>=0);
    assert(l1<=L_MAX);
    assert(l2<=L_MAX);
    assert(l3<=L_MAX);

    assert(m1>=0);
    assert(m2>=0);
    assert(m3>=0);
    assert(m1<=M_COUNT);
    assert(m2<=M_COUNT);
    assert(m3<=M_COUNT);


    double gaunt_coeff = gaunt_table[l1][l2][l3][m1][m2][m3];
    gaunt_coeff *= pow(-1.0, double(orig_m1));

    return gaunt_coeff;

}

}