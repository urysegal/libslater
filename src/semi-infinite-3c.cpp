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


    DECLARE_INDEX_VARIABLE(m_semi_inf)

    indexer_t get_next_sum_from() override { return 1; }

    indexer_t get_next_sum_to() override { return 1; }

    indexer_t get_next_sum_step() override { return 1; }

public:
    Semi_Infinite_Integral_Sum_3(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(
            from_, to_, s, step_) {}


    static complex calculate_expression(Sum_State *state) {
        //compute line 3 and 4 here
        //GAUTAM -> THESE ARE STILL DUMMY VALUES
        auto binomial = bm::binomial_coefficient<double>((2 * state->niu() - state->n_gamma()) / 2.0,
                                                         state->m_semi_inf);

        auto R2S = state->R2 * state->R2 * state->s * (1 - state->s);
        double m_power = pow((R2S * sqrt(state->a()/state->b())) / 2.0, state->m_semi_inf);

        double poch = pochhammer(state->niu() - state->miu(), state->niu() - state->m_semi_inf);

        auto vb = (state->n_x() + state->lambda - state->n_gamma()) / 2.0 + state->sigma + state->m_semi_inf +
                  1.0 / 2.0;
        auto x = sqrt(state->a()/state->b()) * sqrt(R2S + state->v() * state->v());
        double K = bm::cyl_bessel_k(vb, x);
        double denominator = pow(sqrt(R2S + state->v() * state->v()), vb);

        return binomial * m_power * poch * (K / denominator);
    }

};

// line 2 in eqn 56

class Semi_Infinite_Integral_Sum_2 : public Nested_Summation<indexer_t, complex, Semi_Infinite_Integral_Sum_3> {

protected:
    DECLARE_INDEX_VARIABLE(sigma);

    indexer_t get_next_sum_from() override { return 0; }

    indexer_t get_next_sum_to() override { return (2 * STATE->niu() - STATE->n_gamma()) / 2.0; }

    complex expression() override {
        auto s = STATE;
        return calculate_expression(s);
    }

public:

    Semi_Infinite_Integral_Sum_2(int from_, int to_, Summation_State<indexer_t, complex> *s, int step_) : Nested_Summation(
            from_, to_, s, step_) {}

    static complex calculate_expression(Sum_State *s) {
        // compute line 2 expression here
        auto binomial = bm::binomial_coefficient<double>((s->r()), s->sigma);
        auto power = pow((s->v() * s->v() * sqrt(s->a()/s->b())) / 2.0, s->sigma);
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
    Semi_Infinite_Integral_Sum_1(Summation_State<indexer_t, complex> *s) : Nested_Summation(1, 1, s) {}

    //rename it to r = nx - lambda /2

    /// We calculate the expression as public and  static so we can call it directly
    /// from the test suites.
    static complex calculate_expression(Sum_State *state) {
        // compute line 1 expression here
        double power1 = pow(-2.0, state->r());
        double power2 = pow(sqrt(state->a()/state->b()), state->n_x() + state->lambda - state->n_gamma() + 1.0);
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

// functions for r=-1 case
complex glevin(complex s_n, complex w_n, double beta, int n,
               std::vector<complex>& numerator_array,
               std::vector<complex>& denominator_array){
    // The calling program needs to call glevin correctly by starting with n=0
    // glevin computes both the estimates for the numerator and denominator sums for the levin's transformation
    // Starting values for numerator and denominator arrays
    // Adapted from FORTRAN 77 SUBROUTINE GLEVIN in NONLINEAR SEQUENCE TRANSFORMATIONS ..., Weniger et al.

    //Array length check, values populated from 0 to n
    assert(int(numerator_array.size())  >= n+1);
    assert(int(denominator_array.size())>= n+1);

    //Starting values. Computation proceeds backwards from nth index to 0th index.
    numerator_array[n] = s_n / w_n; //7.2-9
    denominator_array[n] = 1.0 / w_n; //7.2-10

    //Inner loop to compute 7.5-3
    if (n>0){
        //factor is just 1 in this case
        //The indexing is not circular, the previous value should be correctly computed in a previous call to glevin for n-1
        numerator_array[n - 1] = numerator_array[n] - numerator_array[n - 1];
        denominator_array[n - 1] = denominator_array[n] - denominator_array[n - 1];
        if (n>1){
            auto bn1 = beta + n - 1;
            auto bn2 = beta + n;
            auto coef = bn1 / bn2;
            for (int j = 2; j<=n; j++){
                auto factor = ( beta + n - j ) * pow(coef,j - 2) / bn2 ;
                //The indexing is not circular, the previous value should be correctly computed in a previous call to glevin for n-1
                numerator_array[n - j] = numerator_array[n - j + 1] - factor * numerator_array[n - j];
                denominator_array[n - j] = denominator_array[n - j + 1] - factor * denominator_array[n - j];
            }
        }
    }

    //Compute estimate of sum= numerator/denominator
    complex sum_est;
//    if (abs(denominator_array[0]) < 1e-10){
//        //Sum is diverging
//        sum_est = 1e60; //RAISE ERROR HERE
//    }
//    else{
    sum_est = numerator_array[0] / denominator_array[0];
//  }
    return sum_est;
}

complex sum_kth_term(Integral_State *state,int p)
{
    auto power = pow(state->beta*state->beta *state->z / 2.0 , p);
    auto pochfrac = 1.0 / pochhammer(state->lambda + 0.5,p+1);
    auto R2S = state->alpha;

    complex sum =0;
    for (int m =0; m <= state->mu; m++){
        auto binomial = bm::binomial_coefficient<double>(state->mu,m);
        auto m_power = pow((R2S * state->z) / 2.0, m);
        double poch = pochhammer(state->nu - state->mu, state->mu - m);
        auto vb = ( state->lambda + state->nu - state->nu+ m + p +  1.0 / 2.0 );
        auto xb = state->z * sqrt(R2S * R2S + state->beta * state->beta);
        double K = bm::cyl_bessel_k(vb, xb);
        double denominator = pow(R2S * R2S + state->beta * state->beta, m/2.0); //outside
        sum = sum + binomial * m_power * poch * K/denominator ;
    }

    //denominator factor pulled out of the inner sum
    double denominator = pow(R2S * R2S + state->beta * state->beta, p/2.0); //outside

    return power * pochfrac * sum/ denominator;
}


complex levin_estimate(Integral_State *state){
    //Initialize variables
    complex s_k = 0;
    complex sum_pre = 0;
    complex sum_cur = 0;
    complex a_k;
    double err_pre = 2.0;
    double err_cur = 1.0;
    int MAX_SUM = 30;

    // Initialize vectors for 0...n values of numerator and denominator
    std::vector<complex> num_array(MAX_SUM+1);
    std::vector<complex> den_array(MAX_SUM+1);

    // Outer Loop for levin's transformation
    // Generate Sequence 7.5-5
    for (int m=0; m<=MAX_SUM; m++){
        //compute ak
        a_k = sum_kth_term(state,m);
        std::cout << m << a_k << "  ";
        //compute sk
        s_k = s_k + a_k;
        //call glevin to get estimate
        sum_cur = glevin(s_k, a_k, 1.0, m, num_array, den_array);
        std::cout << sum_cur << " ";

        err_cur = abs(sum_cur - sum_pre);
        //Do at least 10 iterations before breaking
        if( m > 10 ) {
            //check if a_k is too small, or levin's estimate converged/diverged
            if ( abs(a_k) < 1.e-16 || err_cur < 1e-16 || err_cur > err_pre) {
                break;
            }
        }
        // update previous terms to current terms
        err_pre = err_cur;
        sum_pre = sum_cur;
    }
    std::cout << sum_pre << std::endl;

    return sum_pre;
}

void setup_integral_state(Sum_State *state, Integral_State *i_state) {

    bzero((void *)i_state, sizeof(*i_state));
    i_state->mu = state->miu();
    i_state->nu = state->niu();
    i_state->lambda = state->r();
    i_state->alpha = state->R2*sqrt(vector_length(state->B));
    i_state->beta = state->v();
    i_state->z = sqrt(state->a()/state->b());
}

complex semi_infinite_3c_integral(Sum_State *state)
{
    //Evaluate Integral from top level sum here

    std::cout << std::endl << state->s << std::endl << std::endl;

    complex I;
    if (state->r()==-1){
        //Use Levin Transformation
        Integral_State i_state;
        setup_integral_state(state, &i_state);

        auto sum = levin_estimate(&i_state);

        //denominator factor pulled out of the sums
        auto R2S = i_state.alpha  ;
        auto ordc = (state->lambda + state->miu() - state->niu() +  1.0 / 2.0);
        double denominator = pow(R2S * R2S + state->v() * state->v(), ordc/2.0);

        auto fac1 = 1.0 / pow(state->s*(1-state->s),state->n_gamma()/2.0) ;
        auto fac2 = pow(2.0,state->miu()-1)*pow(i_state.z, state->lambda + state->miu() - state->niu() + 1.0 / 2.0);
        auto fac3 = pow(state->v(), state->lambda);
        I = fac1*fac2*fac3*sum/denominator;
        std::cout << "levin sum " << I << std::endl;

    }
    else{
        // Use formula 58 in :
        // [2]Three-Center Nuclear Attraction_Rv5.pdf
        Semi_Infinite_Integral_Sum_1 top_sum(state);

        I = top_sum.get_value();
    }
    return I;
};

}