#pragma once

#include "libslater.h"
#include "integrators.h"
#include "nested_summation.h"
#include "coordinates.h"

namespace slater {

const std::string analytical_3c_name = "analytical-3c-nuclear-attraction";

typedef int indexer_t; /// For this algorithm, the indices are integers (negative and positive )

/// This is the state of the summation in [1] eqn. 28

struct Sum_State : public Summation_State<indexer_t> {

    /// Problem parameters
    unsigned int n1,n2;
    int l1,l2;
    int m1,m2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;
    center_t A = {};
    center_t B = {};
    center_t C = {};

    // Some precalculated values
    double R2 = 0 ; /// distance between A,B
    center_t R2_point = {}; /// Vector from A to B
    Spherical_Coordinates R2_spherical; /// Spherical coordinate of R2_point above

    /// Nested iteration variables
    indexer_t l1_tag;
    indexer_t m1_tag;
    indexer_t l2_tag;
    indexer_t m2_tag;
    indexer_t l;
    indexer_t lambda; //should change this to lambda!
    indexer_t j;

    /// These are parameters that change every iteration, see [1] eqn. 29
    int n_gamma;
    double v  ; //this was actually redefined in Sum 8 calculate_Ylm, can probably be removed
    int n_x;
    double niu;
    double miu;
    double z;

    ///Parameters for semi-infinite integral, see [1] eqn. 56
    double s;
    int sigma;
    int m_semi_inf;

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
//Semi-Infinite Integral eqn 31 & 56 in [1]
std::complex<double> semi_infinite_integral(Sum_State *state);
void update_dependent_parameters(Sum_State *state);


}