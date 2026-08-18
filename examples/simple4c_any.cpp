#include <iostream>
#include "libslater.h"
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <array>

using namespace std;
using namespace slater;

int main(int argc, char** argv) {
    if (argc != 25) {  // 4 flags for coordinates, each followed by 3 values, and 4 flags for exponents, each followed by 1 value
        cerr << "Usage: " << argv[0] << " -a1 a1 -c1 x1 y1 z1 -a2 a2 -c2 x2 y2 z2 -a3 a3 -c3 x3 y3 z3 -a4 a4 -c4 x4 y4 z4" << endl;
        return 1;
    }

    // Maps to store coordinates and exponents
    unordered_map<string, vector<double>> coordinates_map;
    unordered_map<string, float> exponents_map;

    for (int i = 1; i < argc; i += 6) {  // 6 arguments per set (-a, -c and values)
        // Parse exponents
        string a_flag = argv[i];
        if (a_flag != "-a1" && a_flag != "-a2" && a_flag != "-a3" && a_flag != "-a4") {
            cerr << "Invalid exponent flag: " << a_flag << endl;
            return 1;
        }
        try {
            exponents_map[a_flag] = stof(argv[i + 1]);
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid argument for " << a_flag << ": " << argv[i + 1] << endl;
            return 1;
        } catch (const std::out_of_range& e) {
            cerr << "Out of range argument for " << a_flag << ": " << argv[i + 1] << endl;
            return 1;
        }

        // Parse coordinates
        string c_flag = argv[i + 2];
        if (c_flag != "-c1" && c_flag != "-c2" && c_flag != "-c3" && c_flag != "-c4") {
            cerr << "Invalid coordinate flag: " << c_flag << endl;
            return 1;
        }
        try {
            vector<double> coords = {stod(argv[i + 3]), stod(argv[i + 4]), stod(argv[i + 5])};
            coordinates_map[c_flag] = coords;
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid argument for " << c_flag << ": " << argv[i + 3] << " " << argv[i + 4] << " " << argv[i + 5] << endl;
            return 1;
        } catch (const std::out_of_range& e) {
            cerr << "Out of range argument for " << c_flag << ": " << argv[i + 3] << " " << argv[i + 4] << " " << argv[i + 5] << endl;
            return 1;
        }
    }

    // Convert std::vector<double> to slater::center_t (assuming center_t is std::array<double, 3>)
    auto convert_to_center_t = [](const vector<double>& vec) -> center_t {
        assert(vec.size() == 3);
        return {vec[0], vec[1], vec[2]};
    };

    // Extract the coordinates and exponents for each STO_Basis_Function_Info
    Quantum_Numbers quantum_numbers1 = {1, 0, 0}; // Assuming quantum numbers are not affected by exponents
    Quantum_Numbers quantum_numbers2 = {1, 0, 0};
    Quantum_Numbers quantum_numbers3 = {1, 0, 0};
    Quantum_Numbers quantum_numbers4 = {1, 0, 0};

    STO_Basis_Function_Info fi1(exponents_map["-a1"], quantum_numbers1);
    STO_Basis_Function_Info fi2(exponents_map["-a2"], quantum_numbers2);
    STO_Basis_Function_Info fi3(exponents_map["-a3"], quantum_numbers3);
    STO_Basis_Function_Info fi4(exponents_map["-a4"], quantum_numbers4);

    STO_Basis_Function f1(fi1, convert_to_center_t(coordinates_map["-c1"]));
    STO_Basis_Function f2(fi2, convert_to_center_t(coordinates_map["-c2"]));
    STO_Basis_Function f3(fi3, convert_to_center_t(coordinates_map["-c3"]));
    STO_Basis_Function f4(fi4, convert_to_center_t(coordinates_map["-c4"]));

    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    assert(engine != nullptr);

    STO_Integration_Options parameters;
    engine->init(parameters);

    energy_unit_t result = engine->electron_repulsion({f1, f2, f3, f4});
    std::cout << "4 center result: " << result << std::endl;

    delete engine;
    return 0;
}

