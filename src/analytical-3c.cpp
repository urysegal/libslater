#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>
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
namespace bg = boost::geometry;

using complex = std::complex<double>;
namespace bm = boost::math;

#define STATE (static_cast<Sum_State *> (state))


namespace slater {

    static Analytical_3C_evaluator dummy_3c(analytical_3c_name);

// Sixth line in [1] eqn. 28

    class Sum_8 : public Nested_Summation<indexer_t, complex, Last_Nested_Summation<indexer_t, complex> >
    {

    protected:
        complex expression() override
        {
            auto s = STATE;
            return calculate_expression(s);
        }

        DECLARE_INDEX_VARIABLE(j);

        indexer_t get_next_sum_from() override { return 1; }

        indexer_t get_next_sum_to() override { return 1; }

        indexer_t get_next_sum_step() override { return 1; }


    public:
        Sum_8(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}


        static complex calculate_Ylm(double s, Sum_State *state)
        {
            Quantum_Numbers quantumNumbers({0, state->lambda, state->miu()});


            auto R1 = state->C;

            center_t scaled_R2 = scale_vector(state->R2_point, 1 - s);
            auto v = vector_between(R1, scaled_R2);


            Spherical_Coordinates v_vec_spherical{v};
            auto theta = v_vec_spherical.theta;
            auto phi = v_vec_spherical.phi;

            return eval_spherical_harmonics(quantumNumbers, theta, phi);
        }


        // Seventh line in [1] eqn. 28
        static complex calculate_gaussian_point(const double &s, Sum_State *state) {

            if ( s >  0.23168792 and s < 0.232 ) {
                std::cout << "Stop Here" << std::endl;
            }

            complex result;
            complex power1 = pow(s, state->n2 + state->l2 + state->l1 - state->l1_tag);
            complex power2 = pow(1 - s, state->n1 + state->l1 + state->l2 - state->l2_tag);
            complex ylm = calculate_Ylm(s, state);
            complex prefactor = power1 * power2 * ylm;
            complex semi_inf = calculate_semi_infinite_integral(s, state);
            result = prefactor * semi_inf;
            return result;
        }

        // eqn. 31 in [1]
        static complex calculate_semi_infinite_integral(const double &s, Sum_State *state) {
            state->s = s; //update value of s in State
            state->quad_points.emplace(s);
            complex integral = semi_infinite_3c_integral(state);
            return integral;
        }

        static complex calculate_integral(Sum_State *state) {
            auto f = [&](const double &s) { return calculate_gaussian_point(s, state); };
            state->quad_points.clear();
            complex Q = boost::math::quadrature::gauss<double, 30>::integrate(f, 0, 1);
            std::cout << "-----" << std::endl;
            for ( auto p : state->quad_points) {
                std::cout << p << std::endl;
            }
            std::cout << "-----" << std::endl;
            return Q;
        }

        static complex calculate_expression(Sum_State *s) {
#ifdef SUMMATION_DEBUG
            s->debug_state.pause();
#endif

            complex result = 0;
            auto choose = bm::binomial_coefficient<double>(s->delta_l(), s->j);
            auto numerator = pow(-1, s->j);
            auto x = s->n1 + s->n2 + s->l1 + s->l2 - s->j + 1;
            auto denominator_pow_2 = pow(2, x);
            auto denominator_factorial = bm::factorial<double>(x);
            auto factor = choose * (numerator / (denominator_pow_2 * denominator_factorial));

            complex integral_value = calculate_integral(s);

            result = factor * integral_value;
#ifdef SUMMATION_DEBUG
            s->debug_state.resume();
#endif

            return result;
        }
    };


// Fifth line in [1] eqn. 28

    class Sum_7 : public Nested_Summation<indexer_t, complex, Sum_8> {

    protected:
        virtual complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }


        DECLARE_INDEX_VARIABLE(lambda)

        indexer_t get_next_sum_from() override { return 0; }

        indexer_t get_next_sum_to() override { return STATE->delta_l(); }

        indexer_t get_next_sum_step() override { return 1; }

    public:
        Sum_7(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}


        static complex calculate_expression(Sum_State *s) {
            complex result = 0;
            auto gaunt_part = Gaunt_Coefficient_Engine::get()->calculate({s->l2 - s->l2_tag, s->m2 - s->m2_tag,
                                                                          s->l1 - s->l1_tag, s->m1 - s->m1_tag,
                                                                          s->lambda, s->miu()});

            if (gaunt_part) {
                auto complex_part = pow(-complex{0, 1}, s->lambda);
                result = gaunt_part * complex_part;
            }
            return result;
        }
    };


// Fourth line in [1] eqn. 28

    class Sum_6 : public Nested_Summation<indexer_t, complex, Sum_7> {

    protected:
        virtual complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }


        DECLARE_INDEX_VARIABLE(l)

        indexer_t get_next_sum_from() override {
            return get_l_min(STATE->l1 - STATE->l1_tag,
                             STATE->l2 - STATE->l2_tag,
                             STATE->m1 - STATE->m1_tag,
                             STATE->m2 - STATE->m2_tag);
        }

        indexer_t get_next_sum_to() override { return STATE->l2 - STATE->l2_tag + STATE->l1 - STATE->l1_tag; }

        indexer_t get_next_sum_step() override { return 2; }


    public:
        Sum_6(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

        static indexer_t get_l_min(indexer_t l1, indexer_t l2, indexer_t m1, indexer_t m2) {
            /// Reference [1] eqn. 24
            indexer_t m = std::max(std::abs(l1 - l2), std::abs(m2 - m1));
            int selector = l1 + l2 + m;
            return m + selector % 2;
        }

        static complex calculate_Ylm(Sum_State *s) {
            Quantum_Numbers quantumNumbers({0, s->l, s->m2_tag - s->m1_tag});
            auto theta = s->R2_spherical.theta;
            auto phi = s->R2_spherical.phi;
            return eval_spherical_harmonics(quantumNumbers, theta, phi);
        }

        static complex calculate_expression(Sum_State *s) {
            complex result = 0;

            auto gaunt_part = Gaunt_Coefficient_Engine::get()->calculate({s->l2_tag, s->m2_tag, s->l1_tag, s->m1_tag,
                                                                          s->l, s->m2_tag - s->m1_tag});


            if (gaunt_part) {
                auto radius_part = pow(s->R2, s->l);
                auto ylm_part = calculate_Ylm(s);
                result = gaunt_part * radius_part * ylm_part;
            }
            return result;
        }
    };


// Third line in [1] eqn. 28 , second summation
    class Sum_5 : public Nested_Summation<indexer_t, complex, Sum_6> {

    protected:
        virtual complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }

        DECLARE_INDEX_VARIABLE(m2_tag)

        indexer_t get_next_sum_from() override {
            return Sum_6::get_l_min(STATE->l1_tag, STATE->l2_tag, STATE->m1_tag, STATE->m2_tag);
        }

        indexer_t get_next_sum_to() override { return STATE->l2_tag + STATE->l1_tag; }

        indexer_t get_next_sum_step() override { return 2; }


    public:
        Sum_5(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

        static complex calculate_expression(Sum_State *s) {
            auto coeff = pow(-1.0, s->l2_tag);
            auto g = calculate_gaunt_fraction(s->l2, s->l2_tag, s->m2, s->m2_tag) * pow(-1, s->l2_tag);
            return coeff * g ;
        }

        static complex calculate_gaunt_fraction(indexer_t l, indexer_t l_tag, indexer_t m, indexer_t m_tag) {
            /// Ref [1] eqn. 28 second line
            complex factor = pow(complex(0, 1), l + l_tag);
            //printf("gaunt(%d, %d, %d, %d, %d, %d)\n",l,l_tag,l-l_tag,m,m_tag,m-m_tag);
            complex enumerator = Gaunt_Coefficient_Engine::get()->calculate({l, m, l_tag, m_tag, l - l_tag, m - m_tag});

            complex denominator = bm::double_factorial<double>(2 * l_tag + 1) *
                                  bm::double_factorial<double>(2 * (l - l_tag) + 1);

            return factor * (enumerator / denominator);
        }

    };


// Third line in [1] eqn. 28 , first summation
    class Sum_4 : public Nested_Summation<indexer_t, complex, Sum_5> {

    protected:

        DECLARE_INDEX_VARIABLE(l2_tag)

        indexer_t get_next_sum_from() override { return -1 * STATE->l2_tag; }

        indexer_t get_next_sum_to() override { return STATE->l2_tag; }

    public:
        Sum_4(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
    };


// Second line in [1] eqn. 28 , second summation
    class Sum_3 : public Nested_Summation<indexer_t, complex, Sum_4> {

    protected:
        complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }


        DECLARE_INDEX_VARIABLE(m1_tag)

        indexer_t get_next_sum_from() override { return 0; }

        indexer_t get_next_sum_to() override { return STATE->l2; }

    public:
        Sum_3(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}

        static complex calculate_expression(Sum_State *s) {
            return Sum_5::calculate_gaunt_fraction(s->l1, s->l1_tag, s->m1, s->m1_tag);
        }

    };


// Second line in [1] eqn. 28 , first summation
    class Sum_2 : public Nested_Summation<indexer_t, complex, Sum_3> {

    protected:

        DECLARE_INDEX_VARIABLE(l1_tag)

        indexer_t get_next_sum_from() override { return -1 * STATE->l1_tag; }

        indexer_t get_next_sum_to() override { return STATE->l1_tag; }

    public:
        Sum_2(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(from_, to_, s, step_) {}
    };


// First line in [1] eqn. 28
    class Sum_1 : public Nested_Summation<indexer_t, complex, Sum_2> {

    protected:


        indexer_t get_next_sum_from() override { return 0; }

        indexer_t get_next_sum_to() override { return STATE->l1; }

        complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }

    public:
        Sum_1(Summation_State<indexer_t, complex> *s) : Nested_Summation(1, 1, s) {}

        /// We calculate the expression as public and  static so we can call it directly
        /// from the test suites.
        static complex calculate_expression(Sum_State *s) {

            double numerator =
                    8 * pow(4 * pi, 2) * pow(-1, s->l1 + s->l2)
                    * bm::double_factorial<double>(2 * s->l1 + 1) * bm::double_factorial<double>(2 * s->l2 + 1)
                    * bm::factorial<double>(s->l1 + s->n1 + s->n2 + s->l2 + 1)
                    * pow(s->zeta1, 2 * s->n1 + s->l1 - 1)
                    * pow(s->zeta2, 2 * s->n2 + s->l2 - 1);
            double denominator =
                    bm::factorial<double>(s->n1 + s->l1) * bm::factorial<double>(s->n2 + s->l2);

            return numerator / denominator;
        }

    };


    complex Analytical_3C_evaluator::integrate(const std::vector<STO_Basis_Function> &functions,
                                               const std::vector<center_t> &centers)
   {
        C = centers[0];

        B_functions_representation_of_STO f1(functions[0], functions[0].get_center());
        B_functions_representation_of_STO f2(functions[1], functions[1].get_center());

        energy_unit_t final_result = do_integrate(f1, f2);

        final_result *= functions[0].get_normalization_coefficient() * functions[1].get_normalization_coefficient();

        return final_result;
    }

complex Analytical_3C_evaluator::integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2)
{
        q1 = f1.get_quantum_numbers();
        q2 = f2.get_quantum_numbers();
        zeta1 = f1.get_alpha();
        zeta2 = f2.get_alpha();
        A = f1.get_center();
        B = f2.get_center();
        setup_state();

        auto res = evaluate();

    return res;
}

    complex Analytical_3C_evaluator::evaluate() {
        Sum_1 top_sum(&state);
        return top_sum.get_value();
    }


    void Analytical_3C_evaluator::setup_state() {
        bzero((void *)&state, sizeof(state));
        state.n1 = q1.n;
        state.n2 = q2.n;
        state.l1 = q1.l;
        state.l2 = q2.l;
        state.m1 = q1.m;
        state.m2 = q2.m;
        state.zeta1 = zeta1;
        state.zeta2 = zeta2;
        state.A = A;
        state.B = B;
        state.C = C;

        state.R2 = distance(A, B);
        state.R2_point = vector_between(A, B);
        state.R2_spherical = Spherical_Coordinates(state.R2_point);
    }


    Analytical_3C_evaluator::Analytical_3C_evaluator() : STO_Integrator(3, 1, 1) {}

    void Analytical_3C_evaluator::init(const slater::STO_Integration_Options &params) {
        params.get(Number_of_quadrature_points_Parameter_Name, number_of_quadrature_points);
    }

    STO_Integrator *Analytical_3C_evaluator::clone() const {
        return new Analytical_3C_evaluator();
    }



}