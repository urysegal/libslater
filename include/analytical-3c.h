#pragma once

#include "libslater.h"
#include "integrators.h"
#include "nested_summation.h"

namespace slater {

const std::string analytical_3c_name = "analytical-3c-nuclear-attraction";

typedef int indexer_t; /// For this algorithm, the indices are integers (negative and positive )

/// This is the state of the summation in [1] eqn. 28

struct Sum_State : public Summation_State<indexer_t> {

    /// Problem parameters
    unsigned int n1,n2;
    unsigned int l1,l2;
    int m1,m2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;
    center_t A = {};
    center_t B = {};
    center_t C = {};

    // Some precalculated values
    double R2 = 0 ; /// distance between A,B


    /// Nested iteration variables
    indexer_t l1_tag;
    indexer_t m1_tag;
    indexer_t l2_tag;
    indexer_t m2_tag;
    indexer_t l;
    indexer_t gamma;
    indexer_t j;

    /// These are parameters that change every iteration, see [1] eqn. 29
    int n_gamma ;
    double v  ;
    int u   ;
    int nx   ;
    int delta_l;


    /// This function is called every time the iteration goes forward (at any level of the nesting) to
    /// update some usefull expressions. Can probably be split so we only update the parameters that have changed,
    /// if profiler points to this as a problem.
    void setup_parameters()
    {
        // Follow [1] eqn. 29
        n_gamma = 2*(n1+l1+n2+l2) - (l1_tag + l2_tag) - l + 1;
        v = n1 + n2 + l1 +l2 - l - j + 0.5 ;
        u = (m2 - m2_tag) - ( m1 - m1_tag ) ;
        nx = l1 - l1_tag + l2 - l2_tag ;
        delta_l = (l1_tag + l2_tag - l ) / 2 ;
    }

};


class Analytical_3C_evaluator : public STO_Integrator {

public:

    /// Build an analitical 3-center nuclear attraction integral
    Analytical_3C_evaluator() ;

    /// This constructor is only used for building the factory
    Analytical_3C_evaluator(const std::string &name) : STO_Integrator(name) {}

    /// Release any memory used by the integrator
    virtual ~Analytical_3C_evaluator() = default;

    virtual STO_Integrator *clone() const override;


    /// Initialize the Homeier integrator with a set of options
    /// \param params set of options for the integrator
    virtual void init(const STO_Integration_Options &params) override;

    virtual energy_unit_t integrate(
            const std::vector<STO_Basis_Function> &functions,
            const std::vector<center_t> &centers
    ) override;

    std::complex<double> evaluate();

    void setup_state();

private:

    Sum_State state;

    int number_of_quadrature_points = 30; /// How many quadrature points we should calculate

    Quantum_Numbers q1;
    Quantum_Numbers q2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;
    center_t A = {};
    center_t B = {};
    center_t C = {};

};

}