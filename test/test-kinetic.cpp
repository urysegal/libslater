#include <libslater.h>
#include <math.h>
#include <gaunt.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>
#include "homeier.h"
using namespace slater;

double epsilon = 1e-14;

TEST_CASE( "One Un-normalized Kinetic B function integral", "[kinetic]" ) {

    Homeier_Integrator integrator(integration_types::KINETIC);
    const Quantum_Numbers q1 = {1,0,0};
    const Quantum_Numbers q2 = {1,0,0};

    // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
    center_t c1 = {0, 0, 0};
    center_t c2 = {-2, 0, 0};

    B_function_details f1(1.0, q1, c1);
    B_function_details f2(1.0, q2, c2);
    energy_unit_t S1;
    energy_unit_t S2;
    auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, false);
    std::cout << res<< std::endl;
    CHECK(abs(res.real()  - 7.04871266857357709e-3)< epsilon);
    CHECK(abs(res.imag()  - 0)< epsilon);
    {
        const Quantum_Numbers q1 = {2,0,0};
        const Quantum_Numbers q2 = {1,0,0};

        // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
        center_t c1 = {0, 0, 0};
        center_t c2 = {-2, 0, 0};

        B_function_details f1(1.0, q1, c1);
        B_function_details f2(1.0, q2, c2);
        energy_unit_t S1;
        energy_unit_t S2;
        auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, false);
        std::cout << res<< std::endl;
        CHECK(abs(res.real()  - 4.75788105128716228e-3)< epsilon);
        CHECK(abs(res.imag()  - 0)< epsilon);
    }
    {
        const Quantum_Numbers q1 = {2,0,0};
        const Quantum_Numbers q2 = {2,0,0};

        // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
        center_t c1 = {0, 0, 0};
        center_t c2 = {-2, 0, 0};

        B_function_details f1(1.0, q1, c1);
        B_function_details f2(1.0, q2, c2);
        energy_unit_t S1;
        energy_unit_t S2;
        auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, false);
        std::cout << res<< std::endl;
        CHECK(abs(res.real()  - 3.15429891918667625e-3)< epsilon);
        CHECK(abs(res.imag()  - 0)< epsilon);
    }
}

TEST_CASE( "One Normalized Kinetic B function integral", "[kinetic]" ) {

    Homeier_Integrator integrator(integration_types::KINETIC);
    const Quantum_Numbers q1 = {1,0,0};
    const Quantum_Numbers q2 = {1,0,0};

    // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
    center_t c1 = {0, 0, 0};
    center_t c2 = {-2, 0, 0};

    B_function_details f1(1.0, q1, c1);
    B_function_details f2(1.0, q2, c2);
    energy_unit_t S1;
    energy_unit_t S2;
    auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, true);
    std::cout << res<< std::endl;
    CHECK(abs(res.real()  - 0.11277940269717723)< epsilon);
    CHECK(abs(res.imag()  - 0)< epsilon);
    {
        const Quantum_Numbers q1 = {2,0,0};
        const Quantum_Numbers q2 = {1,0,0};

        // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
        center_t c1 = {0, 0, 0};
        center_t c2 = {-2, 0, 0};

        B_function_details f1(1.0, q1, c1);
        B_function_details f2(1.0, q2, c2);
        energy_unit_t S1;
        energy_unit_t S2;
        auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, true);
        std::cout << res<< std::endl;
        CHECK(abs(res.real()  - 0.11509184026818177)< epsilon);
        CHECK(abs(res.imag()  - 0)< epsilon);
    }
    {
        const Quantum_Numbers q1 = {2,0,0};
        const Quantum_Numbers q2 = {2,0,0};

        // SHIFTED CENTERS FOR THE CASE 2,0,1 and 0,0,1
        center_t c1 = {0, 0, 0};
        center_t c2 = {-2, 0, 0};

        B_function_details f1(1.0, q1, c1);
        B_function_details f2(1.0, q2, c2);
        energy_unit_t S1;
        energy_unit_t S2;
        auto res = integrator.calculate_B_function_kinetic(f1,f2,S1,S2, true);
        std::cout << res<< std::endl;
        CHECK(abs(res.real()  - 0.11535721761596987)< epsilon);
        CHECK(abs(res.imag()  - 0)< epsilon);
    }
}

