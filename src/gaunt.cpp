#include "gaunt.h"



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
    /// Gautam --- please add the real expression
    return args[0] + args[1] * args[3] ;
}

}