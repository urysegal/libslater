#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"

#include <iomanip>

using namespace slater;

int
main()
{
    // -0.254 205 622 595 321×10-1
    center_t c1 = {0, 0, 0};
    center_t c2 = {0,0, 2.5};
    principal_quantum_number_t n1 = 2;
    angular_quantum_number_t l1 = 1;
    magnetic_quantum_number_t m1 = -1;
    sto_exponent_t a1 = 2.4;
    principal_quantum_number_t n2 = 2;
    angular_quantum_number_t l2 = 1;
    magnetic_quantum_number_t m2 = -1;
    sto_exponent_t a2 = 4.1;


//    -0.117 413 789 68 662 828 485 49×10^0
//    center_t c1 = {0, 0, 0};
//    center_t c2 = {0,0, 1.4};
//    principal_quantum_number_t n1 = 1;
//    angular_quantum_number_t l1 = 0;
//    magnetic_quantum_number_t m1 = 0;
//    sto_exponent_t a1 = 10;
//    principal_quantum_number_t n2 = 2;
//    angular_quantum_number_t l2 = 1;
//    magnetic_quantum_number_t m2 = 0;
//    sto_exponent_t a2 = 2;

//    center_t c1 = {1.1000000000000001, 1, 0};
//    center_t c2 = {-1, 1.5, 2.2};
//
//    principal_quantum_number_t n1 = 1;
//    angular_quantum_number_t l1 = 0;
//    sto_exponent_t a1 = 0.70710678118654757;
//
//    principal_quantum_number_t n2 = 3;
//    angular_quantum_number_t l2 = 2;
//    sto_exponent_t a2 = 0.33333334326744080;
    std::cout << "n1,l1,m1,a1: " << n1 << " " << l1 << " " << m1 << " " << a1 << std::endl;
    std::cout << "n2,l2,m2,a2: " << n2 << " " << l2 << " " << m2 << " " << a2 << std::endl;
    std::cout << std::fixed;
    std::cout << std::setprecision(16);
    std::cout << "Overlap: " << std::endl;


    Quantum_Numbers quantum_numbers1 = {n1,l1,m1};
    STO_Basis_Function_Info oxygen_s(  a1 , quantum_numbers1);
    STO_Basis_Function oxygen_1_s(oxygen_s, c1);

    Quantum_Numbers quantum_numbers2 ={n2,l2,m2};
    STO_Basis_Function_Info oxygen_p(a2, quantum_numbers2);
    STO_Basis_Function oxygen_2_p(oxygen_p, c2);
    std::map<slater::integration_types, std::string> engines;
    auto engine = STO_Integration_Engine().create(engines);
    if (engine) {
        STO_Integration_Options options;
        engine->init(options);
        std::vector<STO_Basis_Function> basis_set = {oxygen_1_s, oxygen_2_p};
        auto result = engine->overlap({basis_set[0], basis_set[1]});
        std::cout << result << "  ";
        delete engine;
    }
    std::cout << std::endl;
//    std::cout << "Kinetic: " << std::endl;
//        Quantum_Numbers quantum_numbers1 = {n1,l1,m1};
//        STO_Basis_Function_Info oxygen_s( a1, quantum_numbers1);
//        STO_Basis_Function oxygen_1_s(oxygen_s, c1);
//            Quantum_Numbers quantum_numbers2 ={n2,l2,m2};
//            STO_Basis_Function_Info oxygen_p(a2, quantum_numbers2);
//            STO_Basis_Function_Info hydrogen_s( 1, quantum_numbers1);
//            STO_Basis_Function oxygen_2_p(oxygen_p, c2);
//            std::map<slater::integration_types, std::string> engines;
//            auto engine = STO_Integration_Engine().create(engines);
//            if (engine) {
//                STO_Integration_Options options;
//                engine->init(options);
//                std::vector<STO_Basis_Function> basis_set = {oxygen_1_s, oxygen_2_p};
//                auto result = engine->kinetic({basis_set[0], basis_set[1]});
//                std::cout << result << "  ";
//                delete engine;
//            }
//        std::cout << std::endl;

    return 0;
}
//
// Created by gkluhana on 19/04/23.
//
