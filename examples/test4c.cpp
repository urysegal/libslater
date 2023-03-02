#include <iostream>
#include "libslater.h"
#include <stdlib.h>
#include <iomanip>
#include <assert.h>
#include "../test/three-c-tests.h"

using namespace std;
using namespace slater;

extern struct four_c_tests_t  four_c_tests[]  ;


int
main(int argc, const char *argv[])
{

    if ( argc < 2 ) {
        cerr << "Usage: " << argv[0] << " test-id [to-test-id]" ;
        ::exit(1);
    }

    auto from_test_id = atoi(argv[1]);

    int to_test_id = from_test_id+1;
    if ( argc == 3) {
        to_test_id = atoi(argv[2]);
    }

    auto epsilon = 10e-7;

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    assert(engine != nullptr);
    STO_Integration_Options parameters;
    parameters.set(Use_Normalized_B_Functions_Parameter_Name, true);

    engine->init(parameters);

    for ( auto tid = from_test_id ; tid < to_test_id ; ++ tid ) {

        struct four_c_tests_t &ti = four_c_tests[tid - 1];

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


        double miss = fabs(result.real() - ti.dbar_res);

        if (miss >= epsilon) {
            printf("%d: Fortran: %10.15f    C++:  %10.15f \tmiss: %10.15f \n", tid,
                   ti.dbar_res,
                   result.real(),
                   miss);
        } else {
            printf("%d: Got it (up to epsilon %f) : %10.15f (%10.15f) \n", tid, epsilon, ti.dbar_res,miss);
        }
    }
    delete engine;

}


#include "../test/four-c-tests-data.cc"
