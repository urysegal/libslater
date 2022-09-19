#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"


using namespace slater;

int
main()
{
    STO_Basis_Function_Info oxygen_s(43.5, 0.252, 0);
    STO_Basis_Function_Info oxygen_p(0.5, 1.29872, 1);

    STO_Basis_Function_Info hydrogen_s(3.15, 0.952, 0);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});
    STO_Basis_Function oxygen_2_p(oxygen_p, {0, 0, -0.14142136});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});
    STO_Basis_Function hydrogen_2_s(hydrogen_s, {-0.70710678, 0, 0.56568542});

    std::vector<STO_Basis_Function> basis_set = {oxygen_1_s, oxygen_2_p, hydrogen_1_s, hydrogen_2_s};

    auto engine = STO_Integration_Engine().create("default");
    if (engine) {

        STO_Integration_Options options;
        options.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(options);
        integral_value result = engine->overlap({basis_set[2], basis_set[0]});
        std::cout << result << std::endl;

        delete engine;

    }

	return 0;
}
