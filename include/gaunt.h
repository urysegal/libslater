#pragma once
#include <array>

namespace slater {

/// A singleton engine to calculate/cache Gaunt coefficients

class Gaunt_Coefficient_Engine {

public:

    /// Calculate the value of the Gaunt coefficient at the given arguments
    /// \param args  the six coefficient calculation input,  < l m | l m | l m >
    /// \return Gaunt coefficient
    double calculate(const std::array<int, 6> &args) const;

    static Gaunt_Coefficient_Engine *get();

    static int get_maximal_gaunt_l() ;

private:
    Gaunt_Coefficient_Engine();

    static Gaunt_Coefficient_Engine *instance ;
};


}