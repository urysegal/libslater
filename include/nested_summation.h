#include "libslater.h"

namespace slater {

template<class T, class Next_Summation, int step = 1 , typename indexing_t = int>
class Nested_Summation
{
private:
    Next_Summation inner_summation;
    indexing_t from;
    indexing_t to;

    virtual T scaling_factor() { return dynamic_cast<T>(1); }
    virtual T inner_sum(indexing_t i) = 0;

public:

    Nested_Summation(indexing_t from_, indexing_t to_) : from(from_), to(to_)
    {}

    T get_value()
    {
        T total_sum = dynamic_cast<T>(0);

        for ( indexing_t i = from ; i <= to ; i+= step) {
            total_sum += inner_summation->get_value(i);
        }

        return scaling_factor() * total_sum;
    }

};

}

