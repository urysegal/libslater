#include "gaunt.h"
#include "../precalculate/gaunt-table.h"


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

double Gaunt_Coefficient_Engine::calculate(const std::array<const int, 6> &args) const
{
    auto l1 = args[0];
    auto m1 = args[1] + L_MAX;

    auto l2 = args[2];
    auto m2 = args[3] + L_MAX;

    auto l3 = args[4];
    auto m3 = args[5] + L_MAX;


    return gaunt_table[l1][l2][l3][m1][m2][m3];
}

}