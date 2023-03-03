#include <cmath>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

#include "libslater.h"
using namespace slater;
#include "three-c-tests.h"

extern struct four_c_tests_t  four_c_tests[]  ;

TEST_CASE( "electron repulsion integrals", "[four-center]" )
{

    auto epsilon = 10e-10;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    CHECK(engine != nullptr);

    int misses = 0;

    STO_Integration_Options parameters;
    parameters.set(Quadrature_4C_algorithm_Parameter_Name, std::string("sdbar"));

    for ( int i = 0 ; four_c_tests[i].n1 != 0 ; ++i ) {

        struct four_c_tests_t &ti = four_c_tests[i];

        Quantum_Numbers quantum_numbers1 = {ti.n1, ti.l1, ti.m1};
        Quantum_Numbers quantum_numbers2 = {ti.n2, ti.l2, ti.m2};
        Quantum_Numbers quantum_numbers3 = {ti.n3, ti.l3, ti.m3};
        Quantum_Numbers quantum_numbers4 = {ti.n4, ti.l4, ti.m4};

        STO_Basis_Function_Info fi1(ti.zeta1, quantum_numbers1);
        STO_Basis_Function_Info fi2(ti.zeta2, quantum_numbers2);
        STO_Basis_Function_Info fi3(ti.zeta3, quantum_numbers3);
        STO_Basis_Function_Info fi4(ti.zeta4, quantum_numbers4);

        STO_Basis_Function f1(fi1, {ti.Ax, ti.Ay, ti.Az});
        STO_Basis_Function f2(fi2, {ti.Bx, ti.By, ti.Bz});
        STO_Basis_Function f3(fi3, {ti.Cx, ti.Cy, ti.Cz});
        STO_Basis_Function f4(fi4, {ti.Dx, ti.Dy, ti.Dz});


        energy_unit_t result = engine->electron_repulsion({f1, f2, f3, f4});


        bool miss = true;
        if (fabs(result.real() - ti.dbar_res) >= epsilon) {
            misses++;
        } else {
            miss = false;
        }
        if ( miss ) {
            printf("case %d: Fortran: %f    C++:  %f \n", i,
                   ti.dbar_res,  result.real());
        }
        CHECK(misses==0);

    }

    delete engine;

}

#include "four-c-tests-data.cc"

