#include "nested_summation.h"

namespace slater {

    class Integer_Sum : public Nested_Summation<long, Last_Nested_Summation> {

        virtual long scaling_factor() { return current_value; }
    public:
        Integer_Sum(int from_, int to_) : Nested_Summation(from_, to_)
        {}

    };


    void
    make_sum(long i)
    {
        Integer_Sum ins(0,i);
    }


}