#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"

#include <iomanip>

using namespace slater;

int
main()
{


    Quantum_Numbers quantum_numbers1 = {5,4,-3};
    Quantum_Numbers quantum_numbers2 = {4,3,2};


    STO_Basis_Function_Info oxygen_s( 1, quantum_numbers1);
    STO_Basis_Function_Info oxygen_p(2, quantum_numbers2);

    STO_Basis_Function_Info hydrogen_s( 1, quantum_numbers1);


    STO_Basis_Function oxygen_1_s(oxygen_s, {2, 0, 0});
    STO_Basis_Function oxygen_2_p(oxygen_p, {0, 0, 1});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});
    STO_Basis_Function hydrogen_2_s(hydrogen_s, {-0.70710678, 0, 0.56568542});

    std::vector<STO_Basis_Function> basis_set = {oxygen_1_s, oxygen_2_p, hydrogen_1_s, hydrogen_2_s};

    std::map<slater::integration_types, std::string> engines;
    auto engine = STO_Integration_Engine().create(engines);
    if (engine) {

        STO_Integration_Options options;
        options.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(options);
        auto result = engine->overlap({basis_set[0], basis_set[1]});

        std::cout << std::fixed;
        std::cout << std::setprecision(16);

        std::cout << "Overlap: " << result << std::endl;

        //result = engine->kinetic({basis_set[0], basis_set[1]});
        //std::cout << "Kinetic: " << result << std::endl;


        {

            Quantum_Numbers q1 = {1,0,0};
            Quantum_Numbers q2 = {1,0,0};

            STO_Basis_Function_Info f1_info(1, q1);
            STO_Basis_Function_Info f2_info(1, q2);

            STO_Basis_Function f1(f1_info, {0, 0, 0.0});
            STO_Basis_Function f2(f2_info, {0, 0, 3});


            result = engine->nuclear_attraction({f1,f2}, {0, 0, 2.5});
            std::cout << "Nuclear Attraction: " << result << std::endl;
        }
        delete engine;

    }

	return 0;
}
