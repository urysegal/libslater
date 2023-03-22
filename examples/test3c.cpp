#include <iostream>
#include "libslater.h"
#include <stdlib.h>
#include <iomanip>
#include <assert.h>
#include "../test/three-c-tests.h"

using namespace std;
using namespace slater;

extern struct three_c_tests_t  three_c_tests[]  ;


int
main(int argc, const char *argv[])
{
    if ( argc != 2 ) {
        cerr << "Usage: " << argv[0] << " test-id" ;
        ::exit(1);
    }

    auto epsilon = 10e-7;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    assert(engine != nullptr);


    struct three_c_tests_t &ti = three_c_tests[atoi(argv[1])-1];

    Quantum_Numbers quantum_numbers1 = {ti.n1, ti.l1, ti.m1};
    Quantum_Numbers quantum_numbers2 = {ti.n2, ti.l2, ti.m2};

    STO_Basis_Function_Info fi1( ti.zeta1, quantum_numbers1);
    STO_Basis_Function_Info fi2( ti.zeta2, quantum_numbers2);

    STO_Basis_Function f1(fi1, {ti.Ax, ti.Ay, ti.Az});
    STO_Basis_Function f2(fi2, {ti.Bx, ti.By, ti.Bz});



    STO_Integration_Options parameters;

    engine->init(parameters);

    energy_unit_t result = engine->nuclear_attraction({f1,f2}, {ti.Cx, ti.Cy, ti.Cz});


    double rmiss = fabs(result.real() - ti.result.real) ;
    double imiss = fabs(result.imag() - ti.result.imag) ;
    bool miss = rmiss >= epsilon or imiss >= epsilon;

    if ( miss ) {
            printf("Fortran: %10.15f + %10.15f i   C++:  %10.15f + %10.15f i\trmiss: %10.15f imiss %10.15f\n",
                   ti.result.real, ti.result.imag,
                   result.real(), result.imag(),
                   rmiss, imiss);
    } else {
        printf("Got it (up to epsilon %f) : %10.15f + %10.15f i \n", epsilon, ti.result.real, ti.result.imag);
    }

    delete engine;

}


#include "../test/three-c-tests-data.cc"
