#include <iostream>
#include "libslater.h"
#include <stdlib.h>
#include <assert.h>

using namespace std;
using namespace slater;


int
main()
{
    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    assert(engine != nullptr);
    STO_Integration_Options parameters;

    engine->init(parameters);

    Quantum_Numbers quantum_numbers1 = {1, 0, 0};
    Quantum_Numbers quantum_numbers2 = {1, 0, 0};
    Quantum_Numbers quantum_numbers3 = {1, 0, 0};
    Quantum_Numbers quantum_numbers4 = {1, 0, 0};

    STO_Basis_Function_Info fi1(1, quantum_numbers1);
    STO_Basis_Function_Info fi2(1, quantum_numbers2);
    STO_Basis_Function_Info fi3(1, quantum_numbers3);
    STO_Basis_Function_Info fi4(1, quantum_numbers4);

    STO_Basis_Function f1(fi1, {1, 2, 3});
    STO_Basis_Function f2(fi2, {-2,1, 3});
    STO_Basis_Function f3(fi3, {3,-2, 1});
    STO_Basis_Function f4(fi4, {2, 3, 1});

    energy_unit_t result = engine->electron_repulsion({f1, f2, f3, f4});

   std::cout << "4 center result: " <<  result <<std::endl;

    delete engine;

}
