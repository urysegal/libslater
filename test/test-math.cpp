#include "bfunctions.h"
#include "homeier.h"
#include <math.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

namespace bm = boost::math;
using namespace slater;

TEST_CASE( "Cartesian Coords {3,0,0} to Spherical {3,0,pi/4} ", "[b_func_engine]" ) {
    auto pi = bm::constants::pi<double>();
    center_t r_cart{3,0,0};
    B_function_Engine bFunctionEngine;
    auto r_spher = bFunctionEngine.cartesian_to_spherical(r_cart);
    CHECK(abs(r_spher[0]-3) == 0);
    CHECK(abs(r_spher[2]-pi/2) ==0);
    CHECK(abs(r_spher[1]-0) == 0);
}

TEST_CASE( "shift_first_center_to_origin ", "[homeier]" ) {
    center_t r1{3,1,4};
    center_t r2{4,2,0};
    center_t new_shifted[2];
    shift_first_center_to_origin(r1,r2,new_shifted);
    center_t r1_shifted = new_shifted[0];
    center_t r2_shifted = new_shifted[1];

    CHECK(abs(r1_shifted[0]-0) == 0);
    CHECK(abs(r1_shifted[1]-0) == 0);
    CHECK(abs(r1_shifted[2]-0) == 0);

    CHECK(abs(r2_shifted[0]-1) == 0);
    CHECK(abs(r2_shifted[1]-1) == 0);
    CHECK(abs(r2_shifted[2]+4) == 0);

}

