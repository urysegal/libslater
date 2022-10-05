#include "libslater.h"

namespace slater {

typedef int indexing_t;

class Summation_State
{

};

template<class T, class Next_Summation, int step = 1 >
class Nested_Summation
{
protected:

    indexing_t from;
    indexing_t to;
    Summation_State *state = nullptr;

    indexing_t current_value = 0;

    virtual T scaling_factor() { return 1; }
    virtual indexing_t  get_next_sum_from() { return 0 ;}
    virtual indexing_t  get_next_sum_to() { return -1; }

    virtual void rename_current_value() {}

public:

    Nested_Summation(indexing_t from_, indexing_t to_, Summation_State *state_) : from(from_), to(to_), state(state_)
    {}



    T get_value( )
    {
        T total_sum = 0;


        for ( indexing_t j = from ; j <= to ; j+= step) {

            current_value = j;
            rename_current_value();

            auto scaling = scaling_factor();

            indexing_t next_from = get_next_sum_from();
            indexing_t next_to = get_next_sum_to();
            Next_Summation inner_summation(next_from, next_to, state);

            total_sum +=  scaling * inner_summation.get_value() ;
        }
        return total_sum;
    }

};

template<class T>
class Last_Nested_Summation {
public:
    Last_Nested_Summation(indexing_t from_, indexing_t to_, Summation_State *state_)
    { }

    T get_value() { return 1;}
};


}

