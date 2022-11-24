#include <complex>
#include <boost/math/constants/constants.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>
#include "nested_summation.h"
#include "analytical-3c.h"
#include "slater-utils.h"

namespace bg = boost::geometry;

using complex = std::complex<double>;
namespace bm = boost::math;
#define STATE (static_cast<Sum_State *> (state))

namespace slater {

//Semi-Infinite Integral Sums, eqn 56 in [1]

// line 3 and 4 in eqn 56
    class Semi_Infinite_Integral_Sum_3
            : public Nested_Summation<indexer_t, complex, Last_Nested_Summation<indexer_t, complex> > {

    protected:
        complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }


        indexer_t &get_index_variable() override { return STATE->m_semi_inf; }

        indexer_t get_next_sum_from() override { return 1; }

        indexer_t get_next_sum_to() override { return 1; }

        indexer_t get_next_sum_step() override { return 1; }

    public:
        Semi_Infinite_Integral_Sum_3(int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(
                from_, to_, s, step_) {}


        static complex calculate_expression(Sum_State *state) {
            //compute line 3 and 4 here
            //GAUTAM -> THESE ARE STILL DUMMY VALUES
            auto binomial = bm::binomial_coefficient<double>((2 * state->niu() - state->n_gamma()) / 2.0,
                                                               state->m_semi_inf);

            auto R2S = state->R2 * state->R2 * state->s * (1 - state->s);
            double mpower = pow((R2S * state->z()) / 2.0, state->m_semi_inf);

            double poch = pochhammer(state->niu() - state->miu(), state->niu() - state->m_semi_inf);

            auto vb = (state->n_x() + state->lambda - state->n_gamma()) / 2.0 + state->sigma + state->m_semi_inf +
                      1.0 / 2.0;
            auto x = state->z() * sqrt(R2S + state->v() * state->v());
            double K = bm::cyl_bessel_k(vb, x);
            double denominator = pow(sqrt(R2S + state->v() * state->v()), vb);

            return binomial * mpower * poch * (K / denominator);
        }

    };

// line 2 in eqn 56

    class Semi_Infinite_Integral_Sum_2 : public Nested_Summation<indexer_t, complex, Semi_Infinite_Integral_Sum_3> {

    protected:
        indexer_t &get_index_variable() override { return STATE->sigma; }

        indexer_t get_next_sum_from() override { return 0; }

        indexer_t get_next_sum_to() override { return (2 * STATE->niu() - STATE->n_gamma()) / 2.0; }

        complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }

    public:

        Semi_Infinite_Integral_Sum_2(int from_, int to_, Summation_State<indexer_t> *s, int step_) : Nested_Summation(
                from_, to_, s, step_) {}

        static complex calculate_expression(Sum_State *s) {
            // compute line 2 expression here
            auto binomial = bm::binomial_coefficient<double>((s->r()), s->sigma);
            auto power = pow((s->v() * s->v() * s->z()) / 2.0, s->sigma);
            auto poch = pochhammer((-s->n_x() - s->lambda + 1.0) / 2.0, s->r() - s->sigma);
            double expression =
                    binomial * power * poch;
            return expression;
        }
    };

// line 1 in eqn 56
    class Semi_Infinite_Integral_Sum_1 : public Nested_Summation<indexer_t, complex, Semi_Infinite_Integral_Sum_2> {

    protected:
        indexer_t get_next_sum_from() override { return 0; }

        indexer_t get_next_sum_to() override { return STATE->r(); }

        complex expression() override {
            auto s = STATE;
            return calculate_expression(s);
        }

    public:
        Semi_Infinite_Integral_Sum_1(Summation_State<indexer_t> *s) : Nested_Summation(1, 1, s) {}

        //rename it to r = nx - lambda /2

        /// We calculate the expression as public and  static so we can call it directly
        /// from the test suites.
        static complex calculate_expression(Sum_State *state) {
            // compute line 1 expression here
            double power1 = pow(-2.0, state->r());
            double power2 = pow(sqrt(state->z()), state->n_x() + state->lambda - state->n_gamma() + 1.0);
            double power3 = pow(state->v(), state->lambda + 1.0);
            double numerator =
                    power1 * power2 * power3; // compute line 1 numerator

            double d_power1 = pow(state->R2, 2.0 * state->niu());
            double d_power2 = pow(state->s * (1 - state->s), state->niu() + state->n_gamma() / 2.0);
            double denominator =
                    d_power1 * d_power2; // compute line 1 denominator

            return numerator / denominator;
        }
    };

// functions for =-1 case
complex glevin(complex s_n,complex w_n,double beta, int n, std::vector<complex>& arup,std::vector<complex>& arlo){

    //Insert glevin code here
    arup[n] = s_n / w_n;
    arlo[n] = 1.0 / w_n;
    if (n>0){
        //factor is just 1
        //This is not circular, the previous value should be correctly computed in a previous call to glevin for n=0
        arup[n-1] = arup[n] - arup[n-1];
        arlo[n-1] = arlo[n] - arlo[n-1];
        if (n>1){
            auto bn1 = beta + n-1;
            auto bn2 = beta + n;
            auto coef = bn1/bn2;
            for (int j = 2; j<=n; j++){
                auto factor = (beta+n-j)* pow(coef,j-2) / bn2 ;
                //This is not circular, the previous value should be correctly computed in a previous call to glevin for n=n-1
                arup[n-j] = arup[n-j+1] - factor*arup[n-j];
                arlo[n-j] = arlo[n-j+1] - factor*arlo[n-j];
            }
        }
    }
    //Compute estimate of sum

    complex sum_est;
    if (abs(arlo[0]) < 1e-10){
        //Sum is diverging
        sum_est = 1e60; //RAISE ERROR HERE
    }
    else{
        sum_est = arup[0]/arlo[0];
    }
    return sum_est;
}

complex sum_kth_term(Sum_State *state,int p){
    auto power = pow(state->v()*state->v()*state->z() / 2.0 , p);
    auto pochfrac = 1 / pochhammer(state->lambda + 0.5,p+1);

    complex sum =0;
    for (int m =0; m<=state->miu(); m++){
        auto binomial = bm::binomial_coefficient<double>(state->miu(),m);

        auto R2S = state->R2 * sqrt(state->s * (1 - state->s)) ;
        auto m_power = pow((R2S * state->z()) / 2.0, m);

        double poch = pochhammer(state->niu() - state->miu(), state->miu() - m);

        auto vb = state->lambda + state->miu() - state->niu() + p + m +  1.0 / 2.0;
        auto x = state->z() * sqrt(R2S * R2S + state->v() * state->v());
        double K = bm::cyl_bessel_k(vb, x);

        double denominator = pow(sqrt(R2S * R2S + state->v() * state->v()), vb/2.0);

        sum = sum + binomial * m_power * poch * K/denominator ;
    }

    return power * pochfrac * sum;
}

complex levin_estimate(Sum_State *state){
        //Insert outer loop of levin summation here

        //Initialize variables
        complex s_k = 0;
        complex sum_est_k = 0;
        complex sum_est_kplus1 = 0;
        complex a_k;
        double beta = 1;
        int MAX_SUM = 10;
        std::vector<complex> num_array(MAX_SUM);
        std::vector<complex> den_array(MAX_SUM);
        for (int m=0; m<=MAX_SUM;m++){
            //compute ak
            a_k = sum_kth_term(state,m);
            //compute sk
            s_k = s_k + a_k;
            //call glevin to get estimate
            sum_est_kplus1 = glevin(s_k,a_k,beta,m, num_array, den_array);
            //check convergence
            if (abs(sum_est_kplus1-sum_est_k)< 1e-10 ){
                break;
            }
            sum_est_k = sum_est_kplus1;
        }
        return sum_est_kplus1;
}


}