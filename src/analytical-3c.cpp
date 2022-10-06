#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "libslater.h"
#include "nested_summation.h"

// Comments in this file reference the following works.
// Reference [1]:
// Slevinsky, R.M., Safouhi, H. Compact formulae for three-center nuclear attraction integrals over exponential type functions. J Math Chem 60, 1337–1355 (2022).
// This file implements the method in [1]

auto pi = boost::math::constants::pi<double>();

using complex = std::complex<double>;
namespace bm = boost::math;

#define STATE (static_cast<Sum_State *> (state))

namespace slater {


/// This is the state of the summation in [1] eqn. 28
struct Sum_State : public Summation_State {

    /// Nested iteration variables
    indexing_t l1_tag;
    indexing_t m1_tag;
    indexing_t l2_tag;
    indexing_t m2_tag;
    indexing_t l;
    indexing_t gamma;
    indexing_t j;

    /// These are parameters that change every iteration, see [1] eqn. 29
    int n_gamma ;
    double v  ; // Alexandra - is it really double
    int u   ;
    int nx   ;
    double delta_l;// Alexandra - is it really double ; why the square bracket


    /// Problem parameters
    unsigned int n1,n2;
    unsigned int l1,l2;
    int m1,m2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;

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

// Second line in [1] eqn. 28 , second summation
class Sum_3 : public Nested_Summation<complex , Last_Nested_Summation<complex> > {

protected:
    virtual complex expression() override; // CONTINUE HERE

    virtual indexing_t & get_index_variable() override { return STATE->m1_tag ; }
    virtual indexing_t  get_next_sum_from() override { return 0; }
    virtual indexing_t  get_next_sum_to() override { return STATE->l2; }

public:
    Sum_3( int from_, int to_, Summation_State *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// Second line in [1] eqn. 28 , first summation
class Sum_2 : public Nested_Summation<complex , Sum_3 > {

protected:

    virtual indexing_t & get_index_variable() override { return STATE->l1_tag ; }

    virtual indexing_t  get_next_sum_from() override { return -1 * STATE->l1_tag; }
    virtual indexing_t  get_next_sum_to() override { return STATE->l1_tag; }

public:
    Sum_2( int from_, int to_, Summation_State *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// First line in [1] eqn. 28
class Sum_1 : public Nested_Summation<complex , Sum_2 > {

protected:

    virtual complex expression() override
    {
        auto s = STATE;
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

    virtual indexing_t  get_next_sum_from() override { return 0 ;}
    virtual indexing_t  get_next_sum_to() override { return STATE->l1; }

public:
    Sum_1( Summation_State *s) : Nested_Summation(1, 1, s)
    {}

};


class Analytical_3C_evaluator {

    Sum_State state;

public:
    Analytical_3C_evaluator(const Quantum_Numbers q1_, const Quantum_Numbers q2_, sto_exponent_t zeta1_, sto_exponent_t zeta2_,
                            const center_t &A_, const center_t &B_, const center_t &C_):
        q1(q1_), q2(q2_), zeta1(zeta1_), zeta2(zeta2_), A(A_), B(B_), C(C_)
        {
            setup_convenience_values();
        }

    ~Analytical_3C_evaluator()
    {
    }

    complex evaluate()
    {
        Sum_1 top_sum(&state);
        return top_sum.get_value();
    }

    void setup_convenience_values()
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



private:
    const Quantum_Numbers q1;
    const Quantum_Numbers q2;
    sto_exponent_t zeta1;
    sto_exponent_t zeta2;
    const center_t A;
    const center_t B;
    const center_t C;



};


}
