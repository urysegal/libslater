#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

#include "libslater.h"
using namespace slater;
#include "three-c-tests.h"

extern struct three_c_tests_t  three_c_tests[]  ;

TEST_CASE( "attraction integrals", "[three-center]" )
{

    auto epsilon = 10e-8;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    CHECK(engine != nullptr);

    int imisses = 0;
    int rmisses = 0;
    int count = 0;

    for ( int i = 0 ; three_c_tests[i].n1 != 0 ; ++i , ++count) {

        struct three_c_tests_t &ti = three_c_tests[i];

        Quantum_Numbers quantum_numbers1 = {ti.n1, ti.l1, ti.m1};
        Quantum_Numbers quantum_numbers2 = {ti.n2, ti.l2, ti.m2};

        STO_Basis_Function_Info fi1( ti.zeta1, quantum_numbers1);
        STO_Basis_Function_Info fi2( ti.zeta2, quantum_numbers2);

        STO_Basis_Function f1(fi1, {ti.Ax, ti.Ay, ti.Az});
        STO_Basis_Function f2(fi2, {ti.Bx, ti.By, ti.Bz});



        STO_Integration_Options parameters;
        parameters.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(parameters);

        energy_unit_t result = engine->nuclear_attraction({f1,f2}, {ti.Cx, ti.Cy, ti.Cz});

        bool miss = true;
        if (fabs(result.real() - ti.result.real) >= epsilon) {
            rmisses++;
        } else if (fabs(result.imag() - ti.result.imag) >= epsilon) {
            imisses++;
        } else {
            miss = false;
        }
        if ( miss ) {
            printf("case %d: Fortran: %f + %f i   C++:  %f + %f i\n", i,
                   ti.result.real, ti.result.imag,
                   result.real(), result.imag());
        }


        //CHECK(fabs(result.imag() - ti.result.imag) < epsilon);
        //CHECK(fabs(result.real() - ti.result.real) < epsilon);

    }
    CHECK(rmisses==0);
    CHECK(imisses==0);
    printf("Count %d\n", count);
    delete engine;

}

#include "three-c-tests-data.cc"
