// Minimal STO ERI calculator using libslater.
// Usage: eri  n1 l1 m1 z1 x1 y1 z1   n2 ...   n3 ...   n4 ...   (28 numbers)
// Computes (12|34) = <X1(r1-A) X3(r2-C) | 1/r12 | X2(r1-B) X4(r2-D)>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "libslater.h"

using namespace slater;

int main(int argc, char **argv)
{
    if (argc != 29) {
        std::cerr << "Usage: " << argv[0]
                  << " {n l m zeta x y z} x4   (functions 1,2 = electron 1; 3,4 = electron 2)\n";
        return 1;
    }

    std::array<STO_Basis_Function, 4> fs = {
        STO_Basis_Function({1, {1, 0, 0}}, {0, 0, 0}),
        STO_Basis_Function({1, {1, 0, 0}}, {0, 0, 0}),
        STO_Basis_Function({1, {1, 0, 0}}, {0, 0, 0}),
        STO_Basis_Function({1, {1, 0, 0}}, {0, 0, 0}),
    };
    for (int i = 0; i < 4; ++i) {
        char **a = argv + 1 + 7 * i;
        Quantum_Numbers q = {atoi(a[0]), atoi(a[1]), atoi(a[2])};
        STO_Basis_Function_Info info(atof(a[3]), q);
        fs[i] = STO_Basis_Function(info, {atof(a[4]), atof(a[5]), atof(a[6])});
    }

    STO_Integration_Engine factory;
    std::map<integration_types, std::string> engines;   // empty = defaults
    STO_Integrations *engine = factory.create(engines);

    STO_Integration_Options options;
    options.set(Quadrature_4C_algorithm_Parameter_Name, std::string("sdbar")); // or "WGREP"
    engine->init(options);

    energy_unit_t eri = engine->electron_repulsion(fs);
    std::cout << std::setprecision(12) << eri.real() << "\n";

    delete engine;
    return 0;
}
