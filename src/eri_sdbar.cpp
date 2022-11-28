#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>
#include "eri_sdbar.h"
#include "gaunt.h"


namespace slater {

auto pi = boost::math::constants::pi<double>();
namespace bg = boost::geometry;

using complex = std::complex<double>;
namespace bm = boost::math;

inline auto double_factorial(unsigned i) { return boost::math::double_factorial<double>(i); }
inline auto factorial(unsigned i) { return boost::math::factorial<double>(i); }
inline auto gaunt( int l1, int m1, int l2, int m2, int l3, int m3 )
    { return Gaunt_Coefficient_Engine::get()->calculate({l1, m1, l2, m2, l3, m3}); }

#define STATE (static_cast<SDbar_Sum_State *> (state))


static Electron_Repulsion_SDbar dummy_eri_sdbar(electron_repulsion_sdbar_name);


// Third summation at [4] eqn. 19 , fourth line first summation

class ERI_Sum_4 : public Nested_Summation<indexer_t, complex, Last_Nested_Summation<indexer_t, complex>> {

protected:

    DECLARE_INDEX_VARIABLE(l2_tag)

    indexer_t get_next_sum_from() override { return std::max( -1 * STATE->l2_tag, STATE->m2 - STATE->l2 + STATE->l2_tag ); }

    indexer_t get_next_sum_to() override { return std::min(STATE->l2_tag, STATE->m2 + STATE->l2 - STATE->l2_tag ); }


public:
    ERI_Sum_4(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// Second summation at [4] eqn. 19 third line

class ERI_Sum_3 : public Nested_Summation<indexer_t, complex, ERI_Sum_4> {

protected:

    DECLARE_INDEX_VARIABLE(m1_tag)

    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->l2; }
    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_3(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        complex p = pow( complex(0, 1.0) , s->l1 + s->l1_tag );
        complex numerator = gaunt(s->l1, s->m1, s->l1_tag, s->m1_tag,s->l1 - s->l1_tag, s->m1 - s->m1_tag );
        complex denominator = double_factorial(2*s->l1_tag+1) * double_factorial(2 * ( s->l1 - s->l1_tag ) + 1 );

        return p * (numerator/denominator);
    }
};




// First summation at [4] eqn. 19 third line

class ERI_Sum_2 : public Nested_Summation<indexer_t, complex, ERI_Sum_3> {

protected:

    DECLARE_INDEX_VARIABLE(l1_tag)

    indexer_t get_next_sum_from() override { return std::max( -1 * STATE->l1_tag, STATE->m1 - STATE->l1 + STATE->l1_tag ); }

    indexer_t get_next_sum_to() override { return std::min(STATE->l1_tag, STATE->m1 + STATE->l1 - STATE->l1_tag ); }


public:
    ERI_Sum_2(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};



// initial coefficient at [4] eqn. 19 - up to the first sum (lines 1,2, some of 3)
class ERI_Sum_1 : public Nested_Summation<indexer_t, complex, ERI_Sum_2 > {

protected:


    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->l1; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_1(Summation_State<indexer_t, complex> *s) : Nested_Summation(1, 1, s) {}

    static complex calculate_expression(SDbar_Sum_State *s) {
        // Follow the lines fromr ref. [4] so that it's easy to debug.
        double line1_coeff = 8 * pow(4 * pi, 5) * double_factorial(2*s->l1+1) * double_factorial(2*s->l2+1) ;

        double line1_fraction =
                factorial(s->n1 + s->l1 + s->n2 + s->l2 + 1)
                /
                ( factorial(s->n1 + s->l1) * factorial(s->n2 + s->l2) )
                ;

        double line2_coeff = pow(-1.0, s->l1 + s->l2) * double_factorial(2*s->l3+1) * double_factorial(2*s->l4+1);
        double line2_fraction =
                factorial(s->n3 + s->l3 + s->n4 + s->l4 + 1)
                /
                ( factorial(s->n3 + s->l3) * factorial(s->n4 + s->l4) );
                ;

        double zettas = 1;
        for ( auto i = 0U ; i < 4 ; ++i ) {
            zettas *= pow(s->zeta[i], 2*s->n_as_vec[i] + s->l_as_vec[i] - 1 );
        }

        return line1_coeff * line1_fraction * line2_coeff * line2_fraction * zettas ;
    }

};


Electron_Repulsion_SDbar::Electron_Repulsion_SDbar() : STO_Integrator(3, 1, 1) {}

void Electron_Repulsion_SDbar::init(const slater::STO_Integration_Options &params) {
    params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
}

STO_Integrator *Electron_Repulsion_SDbar::clone() const {
    return new Electron_Repulsion_SDbar();
}


energy_unit_t Electron_Repulsion_SDbar::integrate(const std::vector<STO_Basis_Function> &functions,
                                                  const std::vector<center_t> &centers)
{
    setup_state(functions);
    return evaluate();
}


void Electron_Repulsion_SDbar::setup_state(const std::vector<STO_Basis_Function> &functions)
{

    state.A = functions[0].get_center();
    state.B = functions[1].get_center();
    state.C = functions[2].get_center();
    state.D = functions[3].get_center();
    for ( auto i = 0U ; i < 4 ; ++i ) {
        state.zeta[i] = functions[i].get_exponent();
        state.n_as_vec[i] = functions[i].get_quantum_numbers().n;
        state.l_as_vec[i] = functions[i].get_quantum_numbers().l;
        state.m_as_vec[i] = functions[i].get_quantum_numbers().m;
    }
    state.n1 = state.n_as_vec[0];
    state.n2 = state.n_as_vec[1];
    state.n3 = state.n_as_vec[2];
    state.n4 = state.n_as_vec[3];

    state.l1 = state.l_as_vec[0];
    state.l2 = state.l_as_vec[1];
    state.l3 = state.l_as_vec[2];
    state.l4 = state.l_as_vec[3];

    state.m1 = state.m_as_vec[0];
    state.m2 = state.m_as_vec[1];
    state.m3 = state.m_as_vec[2];
    state.m4 = state.m_as_vec[3];

}

complex Electron_Repulsion_SDbar::evaluate() {
    ERI_Sum_1 top_sum(&state);
    return top_sum.get_value();
}


}