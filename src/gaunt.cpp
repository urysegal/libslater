#include "gaunt.h"



namespace slater {


double Gaunt_Coefficient_Engine::calculate(const std::array<const int, 6> &args) const
{
    /// Gautam --- please add the real expression
    return args[0] + args[1] * args[3] ;
}

}