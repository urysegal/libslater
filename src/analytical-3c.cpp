#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "libslater.h"
#include "nested_summation.h"
#include "gaunt.h"

// Comments in this file reference the following works.
// Reference [1]:
// Slevinsky, R.M., Safouhi, H. Compact formulae for three-center nuclear attraction integrals over exponential type functions. J Math Chem 60, 1337–1355 (2022).
// This file implements the method in [1]

auto pi = boost::math::constants::pi<double>();

using complex = std::complex<double>;
namespace bm = boost::math;

#define STATE (static_cast<Sum_State *> (state))


namespace slater {


typedef int indexer_t; /// For this algorithm, the indices are integers (negative and positive )


/// This is the state of the summation in [1] eqn. 28

struct Sum_State : public Summation_State<indexer_t> {

    /// Problem parameters
    unsigned int n1,n2;
    unsigned int l1,l2;
    int m1,m2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;

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
    double v  ; // Alexandra - is it really double
    int u   ;
    int nx   ;
    double delta_l;// Alexandra - is it really double ; why the square bracket


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

class Analytical_3C_evaluator {

    Sum_State state;

public:
    Analytical_3C_evaluator(const Quantum_Numbers q1_, const Quantum_Numbers q2_, sto_exponent_t zeta1_, sto_exponent_t zeta2_,
                            const center_t &A_, const center_t &B_, const center_t &C_);

    ~Analytical_3C_evaluator() = default;

    complex evaluate();

    void setup_convenience_values();

private:
    const Quantum_Numbers q1;
    const Quantum_Numbers q2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;
    const center_t A;
    const center_t B;
    const center_t C;
};


// Second line in [1] eqn. 28 , second summation
class Sum_5 : public Nested_Summation<indexer_t, complex , Last_Nested_Summation<indexer_t,complex> > {

protected:
    virtual complex expression() override
    {
        auto s = STATE;
        return calculate_expression(s);
    }


    virtual indexer_t & get_index_variable() override { return STATE->m2_tag ; }

    virtual indexer_t  get_next_sum_from() override { return 0; }
    virtual indexer_t  get_next_sum_to() override { return  get_l_min(STATE->l1_tag, STATE->l2_tag, STATE->m1_tag, STATE->m2_tag);}
    virtual indexer_t  get_next_sum_step() override { return 2; }

    indexer_t get_l_min(indexer_t l1, indexer_t l2, indexer_t m1, indexer_t m2)
    {
        /// Reference [1] eqn. 24
        indexer_t m = std::max( std::abs(l1-l2),  std::abs(m2-m1));
        int selector = l1 + l2 + m;
        return m + selector % 2 ;
    }

public:
    Sum_5( int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

    static complex calculate_expression( Sum_State *s )
    {
        s->setup_parameters();
        return calculate_guant_fraction(s->l2, s->l2_tag, s->m2, s->m2_tag) ;
    }

    static complex calculate_guant_fraction(indexer_t l, indexer_t l_tag, indexer_t m, indexer_t  m_tag)
    {
        /// Ref [1] eqn. 28 second line
        complex factor = pow(complex(0,1), l+l_tag) ;
        complex enumerator = Gaunt_Coefficient_Engine::get()->calculate({l,m,l_tag,m_tag, l-l_tag, m-m_tag});

        complex denominator = bm::double_factorial<double>(2*l_tag+1) *  /// WHY SQUARE BRACKET
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
        return Sum_5::calculate_guant_fraction(s->l1, s->l1_tag, s->m1, s->m1_tag) ;
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




Analytical_3C_evaluator::Analytical_3C_evaluator(const Quantum_Numbers q1_, const Quantum_Numbers q2_, sto_exponent_t zeta1_, sto_exponent_t zeta2_,
                                                 const center_t &A_, const center_t &B_, const center_t &C_):
    q1(q1_), q2(q2_), zeta1(zeta1_), zeta2(zeta2_), A(A_), B(B_), C(C_)
{
    setup_convenience_values();
}

complex Analytical_3C_evaluator::evaluate()
{
    Sum_1 top_sum(&state);
    return top_sum.get_value();
}

void Analytical_3C_evaluator::setup_convenience_values()
{
    state.n1=q1.n;
    state.n2=q2.n;
    state.l1=q1.l;
    state.l2=q2.l;
    state.m1=q1.m;
    state.m2=q2.m;
    state.zeta1 = zeta1;
    state.zeta2 = zeta2;
}


}
