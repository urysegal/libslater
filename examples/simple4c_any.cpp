#include <iostream>
#include "libslater.h"
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
using namespace std;
using namespace slater;

int main(int argc, char *argv[]) {
    if (argc != 13) {
        cerr << "Usage: " << argv[0] << " x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4" << endl;
        return 1;
    }

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

    // Parse coordinates from command-line arguments
    vector<double> coordinates; // Changed to double
    for (int i = 1; i < argc; ++i) {
        try {
            coordinates.push_back(stod(argv[i])); // Changed to stod
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid argument: " << argv[i] << endl;
            return 1;
        } catch (const std::out_of_range& e) {
            cerr << "Out of range argument: " << argv[i] << endl;
            return 1;
        }
    }

    STO_Basis_Function f1(fi1, {coordinates[0], coordinates[1], coordinates[2]});
    STO_Basis_Function f2(fi2, {coordinates[3], coordinates[4], coordinates[5]});
    STO_Basis_Function f3(fi3, {coordinates[6], coordinates[7], coordinates[8]});
    STO_Basis_Function f4(fi4, {coordinates[9], coordinates[10], coordinates[11]});

    energy_unit_t result = engine->electron_repulsion({f1, f2, f3, f4});
    std::cout << "4 center result: " <<  result <<std::endl;
    delete engine;
    return 0;
}

