#pragma once

#include "libslater.h"
#define SUMMATION_DEBUG 1

#ifdef SUMMATION_DEBUG
#include "nested_summation_debug.h"
#endif

namespace slater {


#ifdef SUMMATION_DEBUG

#define DECLARE_INDEX_VARIABLE(x) \
    indexer_t &get_index_variable() override { return STATE->x; } \
    const char *get_index_variable_name() override { return #x; }


#else

#define DECLARE_INDEX_VARIABLE(x) \
    indexer_t &get_index_variable() override { return STATE->x; }

#endif

/// Stub class for summation step implementation to keep track of the various summation indices.
template<typename indexing_t, class T>
struct Summation_State
{
    indexing_t _dummy; // to be used if we don't care about the value of the current index.
    unsigned int skipped = 0; /// How many elements were skipped as the coefficient was smaller than the threshold.

#ifdef SUMMATION_DEBUG

    Summation_Debug_State<indexing_t, T> debug_state;

#endif

};

/// One step in a nested summation. Each step is assumed to have some expression, returned by expression(),
/// times another summation, of class "Next_Summation". e.g.  overall value is:
/// sigma(i=from->to) { expression(i)*Next_Summation}
/// You set the from, to and step of the NEXT summation by implementing get_next_sum_from,get_next_sum_to, get_next_sum_step.
/// If inner sum needs to know the value of the indices of the outer sums, subclass Summation_State and implement rename_current_value()
/// so it saved the current value in there according the the nomenclature of your expression.

template<typename indexing_t, class T, class Next_Summation >
class Nested_Summation
{
protected:

    indexing_t from = 1; /// Where this sum starts at
    indexing_t to = 1 ; /// Where this sum stops at. Default from=1, to=1 means to evaluate the expression once.
    indexing_t step = 1; /// How much the index is incremented by
    Summation_State<indexing_t, T> *state = nullptr; /// Point to the state of the whole nested summation so far

    double accuracy_threshold = 0; /// If the coefficient is smaller than that, do not calculate inner sum. override in subclass if desired. Set to 0 to avoid that check.

    virtual indexing_t & get_index_variable() { return state->_dummy; }
    virtual const char *get_index_variable_name() { return "dummy"; }

    /// Calculate the expression for the current itaration (you can fine it in current_index_value). The index values
    /// of the outer summations is accessible in the "state" member (as long as you implement  rename_current_value() and save them
    /// in "state" member... )
    /// \return value for this iteration
    virtual T expression()
    {
        return 1;
    }

    /// The three following funcitons are about the NEXT iteration.


    /// Implement this to control the starting point of the next iteration
    /// \return In the next inner sum, where should the iteration start?
    virtual indexing_t  get_next_sum_from() = 0;

    /// Implement this to control the last of the next iteration.
    /// \return In the next inner sum, what is the last value?
    virtual indexing_t  get_next_sum_to() = 0;

    /// Implement this to control the iteration step size of the next iteration. if you don't, default 1 is used.
    /// \return In the next inner sum, what is the step between iteration?
    virtual indexing_t  get_next_sum_step()
    {
        return 1;
    }


public:

    Nested_Summation(indexing_t from_, indexing_t to_, Summation_State<indexing_t, T> *state_, int step_ = 1) :
        from(from_), to(to_), step(step_), state(state_)
    {}
    virtual ~Nested_Summation() = default;


    /// Calculate the value of the sum, which is the total of expression*inner_sum.
    /// The iteration variable is retured
    /// The expression is retuened by expression() which should be implemented by a subclass, unless
    /// there is just an inner-sum in which case the default returning 1 can be used.
    /// \return total value result of the summation
    T get_value( )
    {
        T total_sum = 0;

        indexing_t & current_index_value = get_index_variable();

#ifdef SUMMATION_DEBUG
        state->debug_state.push_state(typeid(*this).name(), get_index_variable_name(), from, to, step);
#endif
        for ( current_index_value = from ; current_index_value <= to ;current_index_value += step) {

            auto scaling = expression();
            if ( scaling != 0.0 and abs(scaling) > accuracy_threshold  ) {
                Next_Summation inner_summation(get_next_sum_from(), get_next_sum_to(), state, get_next_sum_step());
                auto inner_value = inner_summation.get_value();
                auto new_item = scaling * inner_value;
                total_sum += new_item ;
#ifdef SUMMATION_DEBUG
                state->debug_state.add_item(current_index_value, scaling, inner_value, new_item);
#endif
            } else {
#ifdef SUMMATION_DEBUG
                state->debug_state.add_zero_item(current_index_value, scaling);
#endif
                this->state->skipped++;
            }
        }
#ifdef SUMMATION_DEBUG
        state->debug_state.pop_state(total_sum);
#endif
        return total_sum;
    }

};

template<typename indexing_t, class T>
class Last_Nested_Summation {
public:
    Last_Nested_Summation(indexing_t from_, indexing_t to_, Summation_State<indexing_t, T> *state_, indexing_t step_)
    { }

    T get_value() { return 1;}
};


}

