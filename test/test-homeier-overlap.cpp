#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

#include "libslater.h"
using namespace slater;

struct test_info {
    int n1;
    int l1;
    int m1;
    double alpha;
    double x1;
    double y1;
    double z1;
    int n2;
    int l2;
    int m2;
    double beta;
    double x2;
    double y2;
    double z2;
    double ovlp_real;
    double ovlp_img;
} ;

extern struct test_info  tests[]  ;

TEST_CASE( "overlap integrals", "[homier]" )
{

    auto epsilon = 10e-8;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    CHECK(engine != nullptr);

    for ( int i = 0 ; tests[i].n1 != 0 ; ++i ) {

        struct test_info &ti = tests[i];

        Quantum_Numbers quantum_numbers1 = {ti.n1, ti.l1, ti.m1};
        Quantum_Numbers quantum_numbers2 = {ti.n2, ti.l2, ti.m2};

        STO_Basis_Function_Info fi1( ti.alpha, quantum_numbers1);
        STO_Basis_Function_Info fi2( ti.beta, quantum_numbers2);

        STO_Basis_Function f1(fi1, {ti.x1, ti.y1, ti.z1});
        STO_Basis_Function f2(fi2, {ti.x2, ti.y2, ti.z2});



            STO_Integration_Options parameters;

            engine->init(parameters);

            energy_unit_t result = engine->overlap({f1, f2});
            CHECK(fabs(result.imag() - ti.ovlp_img) < epsilon);
            CHECK(fabs(result.real() - ti.ovlp_real) < epsilon);


    }
        delete engine;

    }

TEST_CASE( "kinetic energy", "[homier]" )
{

    auto epsilon = 10e-8;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    CHECK(engine != nullptr);
    {
        Quantum_Numbers quantum_numbers1 = {1, 0, 0};
        Quantum_Numbers quantum_numbers2 = {1, 0, 0};

        STO_Basis_Function_Info fi1(1, quantum_numbers1);
        STO_Basis_Function_Info fi2(2, quantum_numbers2);

        STO_Basis_Function f1(fi1, {2, 0, 0});
        STO_Basis_Function f2(fi2, {0, 0, 1});


        STO_Integration_Options parameters;

        engine->init(parameters);

        energy_unit_t result = engine->kinetic({f1, f2});
        CHECK(fabs(result.imag() - 0) < epsilon);
        CHECK(fabs(result.real() - 3.16086324190255619e-2) < epsilon);
    }

    {

        Quantum_Numbers quantum_numbers1 = {1,0,0};
        Quantum_Numbers quantum_numbers2 ={1,0,0};


        STO_Basis_Function_Info fi1( 3, quantum_numbers1);
        STO_Basis_Function_Info fi2(2, quantum_numbers2);

        STO_Basis_Function f1(fi1, {2, 0, 1});
        STO_Basis_Function f2(fi2, {0, 0, 1});


        STO_Integration_Options parameters;

        engine->init(parameters);

        energy_unit_t result = engine->kinetic({f1, f2});
        CHECK(fabs(result.imag() - 0) < epsilon);
        CHECK(fabs(result.real() - (-4.07981676840826202e-2)) < epsilon);
    }

    {

        Quantum_Numbers quantum_numbers1 = {1,0,0};
        Quantum_Numbers quantum_numbers2 ={1,0,0};


        STO_Basis_Function_Info fi1( 1, quantum_numbers1);
        STO_Basis_Function_Info fi2(1, quantum_numbers2);

        STO_Basis_Function f1(fi1, {2, 0, 1});
        STO_Basis_Function f2(fi2, {0, 0, 1});


        STO_Integration_Options parameters;

        engine->init(parameters);

        energy_unit_t result = engine->kinetic({f1, f2});
        CHECK(fabs(result.imag() - 0) < epsilon);
        CHECK(fabs(result.real() - 0.11277940269717723) < epsilon);
    }
 
    delete engine;

}

#include "homeier-test-data.h"
