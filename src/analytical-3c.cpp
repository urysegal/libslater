#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "libslater.h"
#include "nested_summation.h"
#include "gaunt.h"
#include "analytical-3c.h"
#include "slater-utils.h"

// Comments in this file reference the following works.
// Reference [1]:
// Slevinsky, R.M., Safouhi, H. Compact formulae for three-center nuclear attraction integrals over exponential type functions. J Math Chem 60, 1337–1355 (2022).
// This file implements the method in [1]

auto pi = boost::math::constants::pi<double>();

using complex = std::complex<double>;
namespace bm = boost::math;

#define STATE (static_cast<Sum_State *> (state))


namespace slater {

static Analytical_3C_evaluator dummy_3c(analytical_3c_name);

// Third line in [1] eqn. 28

class Sum_6 : public Nested_Summation<indexer_t, complex , Last_Nested_Summation<indexer_t,complex> >
{

protected:
    virtual complex expression() override
    {
        auto s = STATE;
        return calculate_expression(s);
    }


    virtual indexer_t & get_index_variable() override { return STATE->l ; }

    virtual indexer_t  get_next_sum_from() override {
        return get_l_min(STATE->l1-STATE->l1_tag,
                         STATE->l2-STATE->l2_tag,
                         STATE->m1-STATE->m1_tag,
                         STATE->m2- STATE->m2_tag);
        }
    virtual indexer_t  get_next_sum_to() override { return STATE->l2 - STATE->l2_tag + STATE->l1 - STATE->l1_tag ;}
    virtual indexer_t  get_next_sum_step() override { return 2; }


public:
    Sum_6( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

    static indexer_t get_l_min(indexer_t l1, indexer_t l2, indexer_t m1, indexer_t m2)
    {
        /// Reference [1] eqn. 24
        indexer_t m = std::max( std::abs(l1-l2),  std::abs(m2-m1));
        int selector = l1 + l2 + m;
        return m + selector % 2 ;
    }

    static complex calculate_Ylm(Sum_State *s)
    {
        Quantum_Numbers quantumNumbers({0, (unsigned  int)s->l, s->m2_tag - s->m1_tag });
        auto theta = s->R2_spherical.theta;
        auto phi = s->R2_spherical.phi;
        return eval_spherical_harmonics(quantumNumbers, theta, phi);
    }

    static complex calculate_expression( Sum_State *s )
    {
        complex result = 0;

        s->setup_parameters();
        auto gaunt_part = Gaunt_Coefficient_Engine::get()->calculate({s->l2_tag, s->m2_tag , s->l1_tag, s->m1_tag,
                                                                     s->l, s->m2_tag-s->m1_tag});


        if ( gaunt_part ) {
            auto radius_part = pow(s->R2, s->l);
            auto ylm_part = calculate_Ylm(s);
            result = gaunt_part * radius_part * ylm_part;
        }
        return result;
    }


};


// Second line in [1] eqn. 28 , second summation
class Sum_5 : public Nested_Summation<indexer_t, complex , Sum_6 >
{

protected:
    virtual complex expression() override
    {
        auto s = STATE;
        return calculate_expression(s);
    }


    virtual indexer_t & get_index_variable() override { return STATE->m2_tag ; }

    virtual indexer_t  get_next_sum_from() override { return Sum_6::get_l_min(STATE->l1_tag, STATE->l2_tag, STATE->m1_tag, STATE->m2_tag); }
    virtual indexer_t  get_next_sum_to() override { return STATE->l2_tag + STATE->l1_tag ;}
    virtual indexer_t  get_next_sum_step() override { return 2; }


public:
    Sum_5( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

    static complex calculate_expression( Sum_State *s )
    {
        s->setup_parameters();
        return calculate_gaunt_fraction(s->l2, s->l2_tag, s->m2, s->m2_tag) * pow(-1, s->l2_tag);
    }

    static complex calculate_gaunt_fraction(indexer_t l, indexer_t l_tag, indexer_t m, indexer_t  m_tag)
    {
        /// Ref [1] eqn. 28 second line
        complex factor = pow(complex(0,1), l+l_tag) ;
        complex enumerator = Gaunt_Coefficient_Engine::get()->calculate({l,m,l_tag,m_tag, l-l_tag, m-m_tag});

        complex denominator = bm::double_factorial<double>(2*l_tag+1) *
                              bm::double_factorial<double>(2*(l - l_tag)+1);

        return factor * ( enumerator / denominator ) ;
    }

};




// Third line in [1] eqn. 28 , first summation
class Sum_4 : public Nested_Summation<indexer_t, complex , Sum_5 > {

protected:

    virtual indexer_t & get_index_variable() override { return STATE->l2_tag ; }

    virtual indexer_t  get_next_sum_from() override { return -1 * STATE->l2_tag; }
    virtual indexer_t  get_next_sum_to() override { return STATE->l2_tag; }

public:
    Sum_4( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// Second line in [1] eqn. 28 , second summation
class Sum_3 : public Nested_Summation<indexer_t, complex , Sum_4 > {

protected:
    virtual complex expression() override
    {
        auto s = STATE;
        return calculate_expression(s);
    }


    virtual indexer_t & get_index_variable() override { return STATE->m1_tag ; }
    virtual indexer_t  get_next_sum_from() override { return 0; }
    virtual indexer_t  get_next_sum_to() override { return STATE->l2; }

public:
    Sum_3( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

    static complex calculate_expression( Sum_State *s )
    {
        s->setup_parameters();
        return Sum_5::calculate_gaunt_fraction(s->l1, s->l1_tag, s->m1, s->m1_tag) ;
    }

};


// Second line in [1] eqn. 28 , first summation
class Sum_2 : public Nested_Summation<indexer_t, complex , Sum_3 > {

protected:

    virtual indexer_t & get_index_variable() override { return STATE->l1_tag ; }

    virtual indexer_t  get_next_sum_from() override { return -1 * STATE->l1_tag; }
    virtual indexer_t  get_next_sum_to() override { return STATE->l1_tag; }

public:
    Sum_2( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// First line in [1] eqn. 28
class Sum_1 : public Nested_Summation<indexer_t, complex , Sum_2 > {

protected:


    virtual indexer_t  get_next_sum_from() override { return 0 ;}
    virtual indexer_t  get_next_sum_to() override { return STATE->l1; }

    virtual complex expression() override
    {
        auto s = STATE;
        return calculate_expression(s);
    }
public:
    Sum_1( Summation_State<indexer_t> *s) : Nested_Summation(1, 1, s)
    {}

    /// We calculate the expression as public and  static so we can call it directly
    /// from the test suites.
    static complex calculate_expression( Sum_State *s )
    {
        s->setup_parameters();

        double enumerator =
            8 * pow ( 4 * pi , 2 ) * pow (-1, s->l1 + s->l2)
            * bm::double_factorial<double>(2*s->l1 + 1) * bm::double_factorial<double>(2*s->l1 + 1)
            * bm::factorial<double>(s->l1 + s->n1 + s->n2 + s->l2  +1 )
            * pow(s->zeta1, 2*s->n1 + s->l1 -1 )
            * pow(s->zeta2, 2*s->n2 + s->l2 -1 );
        double denominator =
            bm::factorial<double>(s->n1 + s->l1) * bm::factorial<double>(s->n2 + s->l2 );

        return enumerator / denominator ;
    }

};






complex Analytical_3C_evaluator::integrate(const std::vector<STO_Basis_Function> &functions,
                                           const std::vector<center_t> &centers)
{
    q1 = functions[0].get_quantum_numbers();
    q2 = functions[1].get_quantum_numbers();
    zeta1 = functions[0].get_exponent();
    zeta2 = functions[1].get_exponent();
    A = functions[0].get_center();
    B = functions[0].get_center();
    C = centers[0];
    setup_state();

    return evaluate();
}


complex Analytical_3C_evaluator::evaluate()
{
    Sum_1 top_sum(&state);
    return top_sum.get_value();
}


void Analytical_3C_evaluator::setup_state()
{
    state.n1=q1.n;
    state.n2=q2.n;
    state.l1=q1.l;
    state.l2=q2.l;
    state.m1=q1.m;
    state.m2=q2.m;
    state.zeta1 = zeta1;
    state.zeta2 = zeta2;
    state.A = A;
    state.B = B;
    state.C = C;

    state.R2 = distance(A, B);
    state.R2_point = vector_between(A, B);
    state.R2_spherical = Spherical_Coordinates(state.R2_point);
}


Analytical_3C_evaluator::Analytical_3C_evaluator() : STO_Integrator(3,1,1)
{}

void Analytical_3C_evaluator::init(const slater::STO_Integration_Options &params)
{
    params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
}

STO_Integrator *Analytical_3C_evaluator::clone() const
{
    return new Analytical_3C_evaluator();
}

}
