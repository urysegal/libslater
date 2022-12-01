#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>
#include "eri_sdbar.h"
#include "gaunt.h"
#include "slater-utils.h"
#include "coordinates.h"

#define STATE (dynamic_cast<SDbar_Sum_State *> (state))
#define miu_from(x) std::max( -1 * STATE->l##x##_tag, STATE->m##x - STATE->l##x + STATE->l##x##_tag )
#define miu_to(x) std::max( STATE->l##x##_tag, STATE->m##x + STATE->l##x - STATE->l##x##_tag )


namespace slater {

auto pi = boost::math::constants::pi<double>();
namespace bg = boost::geometry;

using complex = std::complex<double>;
namespace bm = boost::math;

inline auto double_factorial(unsigned i) { return boost::math::double_factorial<double>(i); }
inline auto factorial(unsigned i) { return boost::math::factorial<double>(i); }
inline auto gaunt( int l1, int m1, int l2, int m2, int l3, int m3 )
    { return Gaunt_Coefficient_Engine::get()->calculate({l1, m1, l2, m2, l3, m3}); }


static indexer_t get_l_min_for_gaunt_summations(indexer_t l1, indexer_t m1, indexer_t l2, indexer_t m2) {
    /// Reference [1] eqn. 24
    indexer_t m = std::max(std::abs(l2 - l1), std::abs(m1 - m2));
    int selector = l1 + l2 + m;
    return m + selector % 2;
}


static Electron_Repulsion_SDbar dummy_eri_sdbar(electron_repulsion_sdbar_name);

// fifteenth summation at [4] eqn. 19 , twelfth line, second summation

class ERI_Sum_15 : public Nested_Summation<indexer_t, complex, Last_Nested_Summation<indexer_t, complex> > {

protected:

    DECLARE_INDEX_VARIABLE(j12)

    indexer_t get_next_sum_from() override { return 0 ;}

    indexer_t get_next_sum_to() override { return STATE->get_delta34(); }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_15(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {

        auto choose1 = bm::binomial_coefficient<double>(s->get_delta12(), s->j12);
        auto choose2 = bm::binomial_coefficient<double>(s->get_delta34(), s->j34);

        auto power = s->get_nu1() + s->get_nu2() + s->l + s->l1_tag + 1 ;

        auto fraction =
                ( pow(-1.0, s->j12 + s->j34))
                /
                ( pow(2.0, power) * factorial(s->get_nu1() + 0.5 + s->l ) * factorial(s->get_nu2() + 0.5 + s->l_tag));

        auto integrals = calculate_integrals(s);

        return choose1 * choose2 * fraction * integrals ;

    }

    static complex calculate_integrals(SDbar_Sum_State *state)
    {
        auto f = [&](const double &s) { return calculate_outer_integral_point(s, state); };
        auto Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);
        return Q;
    }

    static complex calculate_outer_integral_point(const double &s, SDbar_Sum_State *state)
    {
        state->s = s;

        auto coeff =
                ( pow(s, state->n2 + state->l2 + state->l1 )*pow(1.0-s, state->n1 + state->l1 + state->l2) )
                /
                (pow(s, state->l1_tag)  * pow(1.0-s, state->l2_tag) ) ;
        auto integral = calculate_inner_integral(state);
        return coeff * integral ;

    };

    static complex calculate_inner_integral(SDbar_Sum_State *state)
    {

        auto f = [&](const double &t) { return calculate_inner_integral_point(t, state); };
        auto Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);
        return Q;
    }

    static complex calculate_inner_integral_point(const double &t, SDbar_Sum_State *state)
    {
        state->t = t;
        auto coeff =
                ( pow(t, state->n4 + state->l4 + state->l3 )*pow(1.0-t, state->n3 + state->l3 + state->l4) )
                /
                (pow(t, state->l3_tag)  * pow(1.0-t, state->l4_tag) ) ;

        Quantum_Numbers quantumNumbers({0, state->lambda, state->m2 - state->get_miu()});

        Spherical_Coordinates v_vec_spherical{state->v_vec};
        auto theta = v_vec_spherical.theta;
        auto phi = v_vec_spherical.phi;

        auto Y = eval_spherical_harmonics(quantumNumbers, theta, phi  );


        auto semi_infinite = calculate_eri_sdbar_semi_infinite(state);
        return coeff * Y * semi_infinite ;

    };


}

// fourteenth summation at [4] eqn. 19 , twelfth line, first summation

class ERI_Sum_14 : public Nested_Summation<indexer_t, complex, ERI_Sum_15 > {

protected:

    DECLARE_INDEX_VARIABLE(j12)

    indexer_t get_next_sum_from() override { return 0 ;}

    indexer_t get_next_sum_to() override { return STATE->get_delta34(); }


public:
    ERI_Sum_14(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:

};


// thirteen summation at [4] eqn. 19 , eleventh line

class ERI_Sum_13 : public Nested_Summation<indexer_t, complex, ERI_Sum_14 > {

protected:

    DECLARE_INDEX_VARIABLE(lambda)

    indexer_t get_next_sum_from() override { return 0 ;}

    indexer_t get_next_sum_to() override { return STATE->get_delta12(); }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_13(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {

        complex p = pow(complex(0,-1), s->lambda );
        auto g = gaunt(s->l12 , s->get_m21(),
                     s->l34 , s->get_m43(),
                     s->lambda, s->get_miu() );
        return p * g;
    }
};


// twelfth summation at [4] eqn. 19 , tenth line

class ERI_Sum_12 : public Nested_Summation<indexer_t, complex, ERI_Sum_13 > {

protected:

    DECLARE_INDEX_VARIABLE(l34)

    indexer_t get_next_sum_from() override
    {
        return get_l_min_for_gaunt_summations(STATE->l12, STATE->get_m21(), STATE->l34, STATE->get_m43());
    }

    indexer_t get_next_sum_to() override { return STATE->l12 + STATE->l34; }
    indexer_t get_next_sum_step() override { return 2; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_12(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {

        return gaunt(s->l4 - s->l4_tag, s->m4 - s->m4_tag,
                     s->l3 - s->l3_tag, s->m3 - s->m3_tag,
                     s->l34 , s->get_m43() );
    }
};

// eleventh summation at [4] eqn. 19 , ninth line

class ERI_Sum_11 : public Nested_Summation<indexer_t, complex, ERI_Sum_12> {

protected:

    DECLARE_INDEX_VARIABLE(l_tag)

    indexer_t get_next_sum_from() override
    {
        return get_l_min_for_gaunt_summations(STATE->l4 - STATE->l4_tag, STATE->m4 - STATE->m4_tag,
                                              STATE->l3 - STATE->l3_tag, STATE->m3 - STATE->m3_tag);
    }

    indexer_t get_next_sum_to() override { return STATE->l3 - STATE->l3_tag + STATE->l4 -STATE->l4_tag; }
    indexer_t get_next_sum_step() override { return 2; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_11(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        auto g = gaunt(s->l4_tag, s->m4_tag, s->l3_tag, s->m3_tag, s->l_tag , s->m4_tag - s->m3_tag );

        auto r_pow = pow( s->R34 , s->l_tag);

        Quantum_Numbers quantumNumbers({0, s->l_tag, s->m4_tag - s->m3_tag});

        Spherical_Coordinates v_vec_spherical{s->R34_vec};
        auto theta = v_vec_spherical.theta;
        auto phi = v_vec_spherical.phi;

        auto Y = eval_spherical_harmonics(quantumNumbers, theta, phi  );

        return g * r_pow * Y;
    }
};



// tenth summation at [4] eqn. 19 , eighth line

class ERI_Sum_10 : public Nested_Summation<indexer_t, complex, ERI_Sum_11 > {

protected:

    DECLARE_INDEX_VARIABLE(l12)

    indexer_t get_next_sum_from() override
    {
        return get_l_min_for_gaunt_summations(STATE->l4_tag, STATE->m4_tag, STATE->l3_tag, STATE->m3_tag);
    }

    indexer_t get_next_sum_to() override { return STATE->l3_tag + STATE->l4_tag; }
    indexer_t get_next_sum_step() override { return 2; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_10(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {

        return gaunt(s->l2 - s->l2_tag, s->m2 - s->m2_tag,
                     s->l1 - s->l1_tag, s->m1 - s->m1_tag,
                     s->l12 , s->get_m21() );
    }
};

// ninth summation at [4] eqn. 19 , seventh line

class ERI_Sum_9 : public Nested_Summation<indexer_t, complex, ERI_Sum_10> {

protected:

    DECLARE_INDEX_VARIABLE(l)

    indexer_t get_next_sum_from() override
    {
        return get_l_min_for_gaunt_summations(STATE->l2 - STATE->l2_tag, STATE->m2 - STATE->m2_tag,
                                              STATE->l1 - STATE->l1_tag, STATE->m1 - STATE->m1_tag);
    }

    indexer_t get_next_sum_to() override { return STATE->l1 - STATE->l1_tag + STATE->l2 -STATE->l2_tag; }
    indexer_t get_next_sum_step() override { return 2; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_9(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        auto g = gaunt(s->l2_tag, s->m2_tag, s->l1_tag, s->m1_tag, s->l , s->m2_tag - s->m1_tag );

        auto r_pow = pow( s->R21, s->l);

        Quantum_Numbers quantumNumbers({0, s->l, s->m2_tag - s->m1_tag});

        Spherical_Coordinates v_vec_spherical{s->R21_vec};
        auto theta = v_vec_spherical.theta;
        auto phi = v_vec_spherical.phi;

        auto Y = eval_spherical_harmonics(quantumNumbers, theta, phi  );

        return g * r_pow * Y;
    }
};


// eighth summation at [4] eqn. 19 , sixth line second summation

class ERI_Sum_8 : public Nested_Summation<indexer_t, complex, ERI_Sum_9 > {

protected:

    DECLARE_INDEX_VARIABLE(m4_tag)

    indexer_t get_next_sum_from() override
    {
        return get_l_min_for_gaunt_summations(STATE->l2_tag, STATE->m2_tag, STATE->l1_tag, STATE->m1_tag);
    }

    indexer_t get_next_sum_to() override { return STATE->l1_tag + STATE->l2_tag; }
    indexer_t get_next_sum_step() override { return 2; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_8(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        complex p = pow( complex(0, 1.0) , s->l4 + s->l4_tag );
        p *= pow(-1, s->l4_tag);
        complex numerator = gaunt(s->l4, s->m4, s->l4_tag, s->m4_tag,s->l4 - s->l4_tag, s->m4 - s->m4_tag );
        complex denominator = double_factorial(2*s->l4_tag+1) * double_factorial(2 * ( s->l4 - s->l4_tag ) + 1 );

        return p * (numerator/denominator);
    }
};


// seventh summation at [4] eqn. 19 , sixth line first summation

class ERI_Sum_7 : public Nested_Summation<indexer_t, complex, ERI_Sum_8> {

protected:

    DECLARE_INDEX_VARIABLE(l4_tag)

    indexer_t get_next_sum_from() override { return miu_from(4) ; }
    indexer_t get_next_sum_to() override { return miu_to(4) ; }


public:
    ERI_Sum_7(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};




// sixth summation at [4] eqn. 19 , fifth line second summation

class ERI_Sum_6 : public Nested_Summation<indexer_t, complex, ERI_Sum_7> {

protected:

    DECLARE_INDEX_VARIABLE(m3_tag)

    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->l4; }


    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_6(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        complex p = pow( complex(0, 1.0) , s->l3 + s->l3_tag );
        p *= pow(-1, s->l2_tag);
        complex numerator = gaunt(s->l3, s->m3, s->l3_tag, s->m3_tag,s->l3 - s->l3_tag, s->m3 - s->m3_tag );
        complex denominator = double_factorial(2*s->l3_tag+1) * double_factorial(2 * ( s->l3 - s->l3_tag ) + 1 );

        return p * (numerator/denominator);
    }
};


// fifth summation at [4] eqn. 19 , fifth line first summation

class ERI_Sum_5 : public Nested_Summation<indexer_t, complex, ERI_Sum_6> {

protected:

    DECLARE_INDEX_VARIABLE(l3_tag)

    indexer_t get_next_sum_from() override { return miu_from(3) ; }
    indexer_t get_next_sum_to() override { return miu_to(3) ; }


public:
    ERI_Sum_5(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};




// fourth summation at [4] eqn. 19 , fourth line second summation

class ERI_Sum_4 : public Nested_Summation<indexer_t, complex, ERI_Sum_5> {

protected:

    DECLARE_INDEX_VARIABLE(m2_tag)

    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->m2_tag; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_4(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        complex p = pow( complex(0, 1.0) , s->l2 + s->l2_tag );
        p *= pow(-1, s->l2_tag);
        complex numerator = gaunt(s->l2, s->m2, s->l2_tag, s->m2_tag,s->l2 - s->l2_tag, s->m2 - s->m2_tag );
        complex denominator = double_factorial(2*s->l2_tag+1) * double_factorial(2 * ( s->l2 - s->l2_tag ) + 1 );

        return p * (numerator/denominator);
    }
};


// Third summation at [4] eqn. 19 , fourth line first summation

class ERI_Sum_3 : public Nested_Summation<indexer_t, complex, ERI_Sum_4> {

protected:

    DECLARE_INDEX_VARIABLE(l2_tag)

    indexer_t get_next_sum_from() override { return miu_from(2) ; }
    indexer_t get_next_sum_to() override { return miu_to(2) ; }


public:
    ERI_Sum_3(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};


// Second summation at [4] eqn. 19 third line, second summation

class ERI_Sum_2 : public Nested_Summation<indexer_t, complex, ERI_Sum_3> {

protected:

    DECLARE_INDEX_VARIABLE(m1_tag)

    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->l2; }
    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Sum_2(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

private:
    static complex calculate_expression(SDbar_Sum_State *s) {
        complex p = pow( complex(0, 1.0) , s->l1 + s->l1_tag );
        complex numerator = gaunt(s->l1, s->m1, s->l1_tag, s->m1_tag,s->l1 - s->l1_tag, s->m1 - s->m1_tag );
        complex denominator = double_factorial(2*s->l1_tag+1) * double_factorial(2 * ( s->l1 - s->l1_tag ) + 1 );

        return p * (numerator/denominator);
    }
};




// First summation at [4] eqn. 19 third line

class ERI_Sum_1 : public Nested_Summation<indexer_t, complex, ERI_Sum_2> {

protected:

    DECLARE_INDEX_VARIABLE(l1_tag)

    indexer_t get_next_sum_from() override { return miu_from(1) ; }
    indexer_t get_next_sum_to() override { return miu_to(1) ; }


public:
    ERI_Sum_1(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
};



// initial coefficient at [4] eqn. 19 - up to the first sum (lines 1,2, some of 3)
class ERI_Top_Sum : public Nested_Summation<indexer_t, complex, ERI_Sum_1 > {

protected:


    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return STATE->l1; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:
    ERI_Top_Sum(Summation_State<indexer_t, complex> *s) : Nested_Summation(1, 1, s) {}

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



    state.R21_vec = vector_between(state.B, state.A);
    state.R34_vec = vector_between(state.C, state.D);

    state.R21 = vector_length(state.R21_vec);
    state.R34 = vector_length(state.R34_vec);

}

complex Electron_Repulsion_SDbar::evaluate() {
    ERI_Top_Sum top_sum(&state);
    return top_sum.get_value();
}


}