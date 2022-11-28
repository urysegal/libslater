#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"

#include <iomanip>

using namespace slater;

int
main()
{
    std::map<slater::integration_types, std::string> engines;
    auto engine = STO_Integration_Engine().create(engines);

    if (engine) {

        STO_Integration_Options options;

        engine->init(options);
        Quantum_Numbers q1 = {4,2,2};
        Quantum_Numbers q2 = {4,2,2};
        Quantum_Numbers q3 = {4,2,2};
        Quantum_Numbers q4 = {4,2,2};

        STO_Basis_Function_Info f1_info(1.5, q1);
        STO_Basis_Function_Info f2_info(0.5, q2);
        STO_Basis_Function_Info f3_info(2.5, q3);
        STO_Basis_Function_Info f4_info(3.5, q4);

        STO_Basis_Function f1(f1_info, {0, 0, 0.0});
        STO_Basis_Function f2(f2_info, {0, 0, 2.5});
        STO_Basis_Function f3(f3_info, {0, 0, 4.5});
        STO_Basis_Function f4(f4_info, {0, 0, 1.5});


        auto result = engine->electron_repulsion({f1, f2, f3, f4});
        std::cout << "Nuclear Attraction: " << result << std::endl;

        delete engine;

    }

    return 0;
}
