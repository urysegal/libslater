#include "libslater.h"

namespace slater {

typedef int indexing_t;


/// Stub class for summation step implementation to keep track of the various summation indices.
class Summation_State
{
};

/// One step in a nested summation. Each step is assumed to have some expression, returned by scaling_factor(),
/// times another summation, of class "Next_Summation". e.g.  overall value is:
/// sigma(i=from->to) { scaling_factor(i)*Next_Summation}
/// You set the from, to and step of the NEXT summation by implementing get_next_sum_from,get_next_sum_to, get_next_sum_step.
/// If inner sum needs to know the value of the indices of the outer sums, subclass Summation_State and implement rename_current_value()
/// so it saved the current value in there according the the nomenclature of your expression.

template<class T, class Next_Summation >
class Nested_Summation
{
protected:

    indexing_t from = 0; /// Where this sum starts at
    indexing_t to = -1 ; /// Where this sum stops at. Default from=-1, to=0 means no work will take place
    indexing_t step = 1; /// How much the index is incremented by
    Summation_State *state = nullptr; /// Point to the state of the whole nested summation so far

    indexing_t current_index_value = 0; /// In this summation, current index value

    /// Calculate the expression for the current itaration (you can fine it in current_index_value). The index values
    /// of the outer summations is accessible in the "state" member (as long as you implement  rename_current_value() and save them
    /// in "state" member... )
    /// \return value for this iteration
    virtual T scaling_factor()
    {
        return 1;
    }

    ///


    /// Implement this to control the starting point of the next iteration
    /// \return In the next inner sum, where should the iteration start?
    virtual indexing_t  get_next_sum_from()
    {
        return 0 ;
    }

    /// Implement this to control the last of the next iteration.
    /// \return In the next inner sum, what is the last value?
    virtual indexing_t  get_next_sum_to() {
        return -1;
    }

    /// Implement this to control the iteration step size of the next iteration. if you don't, default 1 is used.
    /// \return In the next inner sum, what is the step between iteration?
    virtual indexing_t  get_next_sum_step()
    {
        return 1;
    }

    virtual void rename_current_value()
    {

    }

public:

    Nested_Summation(indexing_t from_, indexing_t to_, Summation_State *state_, int step_ = 1) :
        from(from_), to(to_), step(step_), state(state_)
    {}



    T get_value( )
    {
        T total_sum = 0;


        for ( current_index_value = from ; current_index_value <= to ;current_index_value += step) {

            rename_current_value();

            auto scaling = scaling_factor();
            Next_Summation inner_summation(get_next_sum_from(), get_next_sum_to(), state, get_next_sum_step());

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

