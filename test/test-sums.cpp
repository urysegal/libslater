#include "nested_summation.h"

#include <math.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

using namespace slater;

struct sum_state : public Summation_State {
    int i=0;
    int j=0;
};


class Integer_Sum : public Nested_Summation<long, Last_Nested_Summation<long> > {

protected:

    virtual long expression() { return static_cast<sum_state *> (state)->j; }

    virtual indexing_t & get_index_variable() { return static_cast<sum_state *> (state)->j; }

public:
    Integer_Sum(int from_, int to_, Summation_State *s, int step_ = 1) : Nested_Summation(from_, to_, s, step_)
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

    class Integer_Sum_Sum : public Nested_Summation<long, Integer_Sum> {

    protected:

        virtual long expression() { return static_cast<sum_state *> (state)->i; }

        virtual indexing_t & get_index_variable() { return static_cast<sum_state *> (state)->i; }

        virtual indexing_t  get_next_sum_from() { return 0 ;}
        virtual indexing_t  get_next_sum_to() { return static_cast<sum_state *> (state)->i; }

    public:
        Integer_Sum_Sum(int from_, int to_, Summation_State *s) : Nested_Summation(from_, to_, s)
        {}

    };


    sum_state s;
    Integer_Sum_Sum ins(0,5,&s);

    CHECK(ins.get_value() == 140 );
}