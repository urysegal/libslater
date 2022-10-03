#pragma once
#include <array>

namespace slater {

class Gaunt_Coefficient_Engine {

public:

    /// Calculate the value of the Gaunt coefficient at the given arguments
    /// \param args  the six coefficient calculation input,  < l m | l m | l m >
    /// \return Gaunt coefficient
    double calculate(const std::array<const int, 6> &args) const;

private:

};


}