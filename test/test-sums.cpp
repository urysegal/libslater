#include "nested_summation.h"

#include <math.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>
#include "analytical-3c.h"

using namespace slater;
using complex = std::complex<double>;

struct sum_state : public Summation_State<int, long> {
    int i=0;
    int j=0;
};


class Integer_Sum : public Nested_Summation<int, long, Last_Nested_Summation<int, long> > {

protected:

    virtual long expression() { return static_cast<sum_state *> (state)->j; }

    virtual int & get_index_variable() { return static_cast<sum_state *> (state)->j; }

    virtual int  get_next_sum_from() { return 1; }

    virtual int  get_next_sum_to() { return 1; }

public:
    Integer_Sum(int from_, int to_, Summation_State<int, long> *s, int step_ = 1) : Nested_Summation(from_, to_, s, step_)
    {}

};

TEST_CASE( "simple sum", "[sums]" )
{
    sum_state s;
    Integer_Sum ins(0,5,&s);

    CHECK(ins.get_value() == 1+2+3+4+5 );
}

TEST_CASE( "simple sum with step", "[sums]" ) {

    sum_state s;
    Integer_Sum ins(3,15,&s,2);

    CHECK(ins.get_value() == 3+5+7+9+11+13+15 );
}



TEST_CASE( "two simple sums", "[sums]" ) {

    class Integer_Sum_Sum : public Nested_Summation<int, long, Integer_Sum> {

    protected:

        virtual long expression() { return static_cast<sum_state *> (state)->i; }

        virtual int & get_index_variable() { return static_cast<sum_state *> (state)->i; }

        virtual int  get_next_sum_from() { return 0 ;}
        virtual int  get_next_sum_to() { return static_cast<sum_state *> (state)->i; }

    public:
        Integer_Sum_Sum(int from_, int to_, Summation_State<int, long> *s) : Nested_Summation(from_, to_, s)
        {}

    };


    sum_state s;
    Integer_Sum_Sum ins(0,5,&s);

    CHECK(ins.get_value() == 140 );
}

TEST_CASE( "alternating series sum", "[sums]" ) {
//Initialize variables
    complex s_k = 0;
    complex sum_est_k = 0;
    complex sum_est_kplus1 = 0;
    complex a_k;
    int MAX_SUM = 100;
    std::vector<complex> num_array(MAX_SUM+1);
    std::vector<complex> den_array(MAX_SUM+1);
    // Outer Loop for levin's transformation
    // Generate Sequence for alternating harmonic series
    for (int m=0; m<=MAX_SUM;m++){
        //compute ak
        a_k = pow(-1.0,m+1) / double(m+1);
        //compute sk
        s_k = s_k + a_k;
        //call glevin to get estimate
        sum_est_kplus1 = glevin(s_k,a_k, 1.0,m, num_array, den_array);
        //check convergence
        if (abs(sum_est_kplus1-sum_est_k) < 1e-17 ){
            break;
        }
        sum_est_k = sum_est_kplus1;
    }
    CHECK(abs(sum_est_kplus1 - -0.6931471805599453) < 1e-14 );

}