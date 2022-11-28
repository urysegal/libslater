#include "nested_summation.h"

namespace slater {


#ifdef SUMMATION_DEBUG

constexpr const int tab_size = 2;
void print_tabs(int t)
{
    for ( auto i = 0 ; i < t*tab_size ; ++i ) {
        std::cerr << " ";
    }
}


#endif


}