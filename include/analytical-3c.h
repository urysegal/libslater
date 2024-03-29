#pragma once

#include "libslater.h"
#include "integrators.h"
#include "nested_summation.h"
#include "coordinates.h"
#include "slater-utils.h"
#include "bfunctions.h"

using complex = std::complex<double>;

namespace slater {

const std::string analytical_3c_name = "analytical-3c-nuclear-attraction";

typedef int indexer_t; /// For this algorithm, the indices are integers (negative and positive )

//
struct Sum_State;
struct Integral_State;

/// Semi-Infinite Integral eqn 31 & 56 in [1] ( see the implementation code file and Readme )
/// \param state state of the nested summation
/// \return value of the semi infinite integral in [1]
std::complex<double> semi_infinite_3c_integral(Sum_State *state);

/// Glevin function adapted from FORTRAN 77 SUBROUTINE GLEVIN in NONLINEAR SEQUENCE TRANSFORMATIONS ..., Weniger et al.
/// \param s_n
/// \param w_n
/// \param beta
/// \param n
/// \param numerator_array
/// \param denominator_array
/// \return L_m^0
complex glevin(complex s_n, complex w_n, double beta, int n,
                   std::vector<complex>& numerator_array,
                   std::vector<complex>& denominator_array);


/// This is the state of the summation in [1] eqn. 28
struct Sum_State : public Summation_State<indexer_t, std::complex<double> > {

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

    /// Nested iteration variables, see [1] eqn. 28
    indexer_t l1_tag;
    indexer_t m1_tag;
    indexer_t l2_tag;
    indexer_t m2_tag;
    indexer_t l;
    indexer_t lambda;
    indexer_t j;


    /// These are parameters that change every iteration, see [1] eqn. 29
    int n_gamma() const{
        return 2*(n1+l1+n2+l2)-(l1_tag+l2_tag)-l+1;
    }
    double niu() const {
        return n1 + n2 + l1 + l2 -l - j + 1.0/2.0;
    }
    int miu() const{
        return (m2-m2_tag) - (m1 - m1_tag);
    }
    double n_x() const{
      return   l1 - l1_tag + l2 - l2_tag;
    };
    double delta_l() const{
        return (l1_tag + l2_tag -l)/2.0;
    }
    double v() const{
        auto R1 = C;
        center_t scaled_R2 =  scale_vector( R2_point, 1-s);
        auto v_vec = vector_between(R1,scaled_R2 );
        return vector_length(v_vec);
    }

    ///Integration variable for semi-infinite integral, see[1] eqn. 28, lines 7-8
    double s;


    /// These are parameters for semi-infinite integral that change every iteration, see [1] eqn. 55

    double a() const{
        return  (1.0-s)*zeta1*zeta1 + s*zeta2*zeta2 ;
    }
    double b() const{
        return s*(1.0-s);
    }
    double r() const{
        return (n_x()-lambda)/2.0 -1;
    }
public:
    ~Sum_State() override = default;

};

struct Integral_State : public Summation_State<indexer_t, std::complex<double> > {
    double mu;
    double nu;
    indexer_t lambda;
    double r;
    double alpha;
    double beta;
    double z;

    int inu; // integer rep of niu

    ///Nested iteration variables for semi-infinite integral, see [1] eqn. 56
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
    ~Analytical_3C_evaluator() override = default;

    STO_Integrator *clone() const override;


    /// Initialize the Homeier integrator with a set of options
    /// \param params set of options for the integrator
    void init(const STO_Integration_Options &params) override;

    energy_unit_t integrate(
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
    sto_exponent_t zeta1 = 0;
    sto_exponent_t zeta2 = 0;
    center_t A = {};
    center_t B = {};
    center_t C = {};

    [[nodiscard]] complex integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2) override;


};


}