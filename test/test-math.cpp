#include "bfunctions.h"
#include "homeier.h"
#include <math.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include "coordinates.h"

namespace bm = boost::math;
using namespace slater;

TEST_CASE( "Cartesian Coords {3,0,0} to Spherical {3,0,pi/4} ", "[b_func_engine]" ) {
    auto pi = bm::constants::pi<double>();
    center_t r_cart{3,0,0};

    Spherical_Coordinates r_spher(r_cart);
    CHECK(abs(r_spher.radius-3) == 0);
    CHECK(abs(r_spher.phi-pi/2) ==0);
    CHECK(abs(r_spher.theta-0) == 0);
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

TEST_CASE( "Evaluate Spherical Harmonics ", "[b_func_engine]" ) {
    auto pi = bm::constants::pi<double>();
    Quantum_Numbers q1 = {4,2,1};
    Quantum_Numbers q2 = {4,2,-1};
    double theta = pi/4;
    double phi = 0;
    B_function_Engine b_func_engine;
    auto Y1 = b_func_engine.eval_spherical_harmonics(q1,theta,phi);
    auto Y2 = b_func_engine.eval_spherical_harmonics(q2,theta,phi);
    CHECK(abs(Y1 - std::complex<double>(-0.3862742020231,0)) < 0.001) ;
    CHECK(abs(Y2 - std::complex<double>(0.3862742020231,0)) < 0.001) ;

}

TEST_CASE( "Evaluate B Functions ", "[b_func_engine]" ) {
    ///TEST MAY NEED FIXING
    Quantum_Numbers q1 = {4,2,1};
    double alpha = 1;
    center_t r{1,2,1} ;
    B_function_Engine b_func_engine;
    auto B = b_func_engine.calculate(q1,alpha,r);
    CHECK(abs(B - std::complex<double>(0.000324863191,0.0001452832357)) < 0.001) ;

}
