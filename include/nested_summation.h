#include "libslater.h"

namespace slater {



template<class T, class Next_Summation, int step = 1 , typename indexing_t = int>
class Nested_Summation
{
protected:

    indexing_t from;
    indexing_t to;
    indexing_t current_value = 0;

    virtual T scaling_factor(indexing_t) { return 1; }
    virtual indexing_t  get_next_sum_from() { return 0 ;}
    virtual indexing_t  get_next_sum_to() { return -1; }
    virtual indexing_t  get_next_sum_step() { return 1; }

    virtual void rename_current_value() {}

public:

    Nested_Summation(indexing_t from_, indexing_t to_) : from(from_), to(to_)
    {}

    T get_value(indexing_t i)
    {
        T total_sum = 0;

        current_value = i;
        rename_current_value();

        indexing_t next_from = get_next_sum_from();
        indexing_t next_to = get_next_sum_to();
        indexing_t next_sum_step = get_next_sum_step();

        Next_Summation inner_summation(next_from, next_to);
        for ( indexing_t j = next_from ; j <= next_to ; j+= next_sum_step) {
            auto scaling = scaling_factor(j);
            total_sum +=  scaling * inner_summation->get_value(j) ;
        }
        return total_sum;
    }

};

class Last_Nested_Summation {
public:

};


}

