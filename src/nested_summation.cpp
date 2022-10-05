#include "nested_summation.h"

namespace slater {


    struct sum_state : public Summation_State {
        int i ;
        int j;
    };

    class Integer_Sum : public Nested_Summation<long, Last_Nested_Summation<long> > {

    protected:

        virtual long scaling_factor() { return static_cast<sum_state *> (state)->j; }

        virtual void rename_current_value() { static_cast<sum_state *> (state)->j = current_value ;}

    public:
        Integer_Sum(int from_, int to_, Summation_State *s) : Nested_Summation(from_, to_, s)
        {}

    };

class Integer_Sum_Sum : public Nested_Summation<long, Integer_Sum> {

protected:

    virtual long scaling_factor() { return static_cast<sum_state *> (state)->i; }

    virtual void rename_current_value() { static_cast<sum_state *> (state)->i = current_value ;}

    virtual indexing_t  get_next_sum_from() { return 0 ;}
    virtual indexing_t  get_next_sum_to() { return static_cast<sum_state *> (state)->i; }

public:
    Integer_Sum_Sum(int from_, int to_, Summation_State *s) : Nested_Summation(from_, to_, s)
    {}

};



    long
    make_sum(long i)
    {
        sum_state s;
        Integer_Sum_Sum ins(0,i,&s);
        return ins.get_value();
    }


}