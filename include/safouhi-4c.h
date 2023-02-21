#pragma once
#include <string>
#include "libslater.h"
#include "integrators.h"

namespace slater {

const std::string safouhi_4c_name = "safouhi-4c-electron-repulsion";

class Safouhi_4C_evaluator : public STO_Integrator {

public:

    /// Build an analitical 3-center nuclear attraction integral
    Safouhi_4C_evaluator() ;

    /// This constructor is only used for building the factory
    Safouhi_4C_evaluator(const std::string &name) : STO_Integrator(name) {}

    /// Release any memory used by the integrator
    ~Safouhi_4C_evaluator() override = default;

    STO_Integrator *clone() const override;


    /// Initialize the 4C integrator with a set of options
    /// \param params set of options for the integrator
    void init(const STO_Integration_Options &params) override;

    energy_unit_t integrate(
        const std::vector<STO_Basis_Function> &functions,
        const std::vector<center_t> &centers
    ) override;


private:

    Quantum_Numbers q1;
    Quantum_Numbers q2;
    Quantum_Numbers q3;
    Quantum_Numbers q4;

    sto_exponent_t zeta1 = 0;
    sto_exponent_t zeta2 = 0;
    sto_exponent_t zeta3 = 0;
    sto_exponent_t zeta4 = 0;

    center_t A = {};
    center_t B = {};
    center_t C = {};
    center_t D = {};

    [[nodiscard]] std::complex<double> integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2) override;


};

}