#include "bfunctions.h"
#include "homeier.h"
#include <math.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include "coordinates.h"
#include "slater-utils.h"

namespace bm = boost::math;
using namespace slater;

TEST_CASE("3D points distance","[utils]")
{
    center_t p0{0, 0, 0};
    center_t p1{1, 1, 1};
    center_t p2{0, 1, 0};

    CHECK(distance(p2,p0) == 1);
    CHECK(abs(distance(p2,p1) - 1.4142135624) < 10E-9);
    CHECK(abs(distance(p0,p1) - 1.73205080756) < 10E-9);
}


TEST_CASE("Angle between points","[utils]")
{
    auto pi = bm::constants::pi<double>();

    center_t p0{0, 0, 0};
    center_t p1{1, 1, 1};
    center_t p2{0, 1, 0};

    center_t connector = vector_between(p0, p1);
    CHECK(p1 == connector);

    center_t connector2 = vector_between(p1, p2);
    CHECK(center_t{ -1,0,-1}  == connector2);

    Spherical_Coordinates sc(connector);
    CHECK(abs(sc.theta - 0.9553166181) < 10E-9 );
    CHECK(sc.phi == pi/4);
}

TEST_CASE( "Cartesian Coords  to Spherical", "[utils]" ) {
    auto pi = bm::constants::pi<double>();

    {
    center_t r_cart{3, 0, 0};
    Spherical_Coordinates r_spher(r_cart);
    CHECK(abs(r_spher.radius - 3) == 0);
    CHECK(abs(r_spher.theta - pi / 2) == 0);
    CHECK(abs(r_spher.phi - 0) == 0);
    }

    {
        center_t r_cart2{0, 0, 1};
        Spherical_Coordinates r_spher(r_cart2);
        CHECK(r_spher.radius == 1);
        CHECK(r_spher.phi == 0);
        CHECK(r_spher.theta == 0);
    }
    {
        center_t r_cart2{1, 0, 0};
        Spherical_Coordinates r_spher(r_cart2);
        CHECK(r_spher.radius == 1);
        CHECK(r_spher.phi == 0);
        CHECK(r_spher.theta == pi/2);
    }
    {
        center_t r_cart2{-1, 0, 0};
        Spherical_Coordinates r_spher(r_cart2);
        CHECK(r_spher.radius == 1);
        CHECK(r_spher.phi == pi);
        CHECK(r_spher.theta == pi/2);
    }
    {
        center_t r_cart2{0, -1, 0};
        Spherical_Coordinates r_spher(r_cart2);
        CHECK(r_spher.radius == 1);
        CHECK(r_spher.phi == -pi/2);
        CHECK(r_spher.theta == pi/2);
    }
    {
        center_t r_cart2{0, 0, -1};
        Spherical_Coordinates r_spher(r_cart2);
        CHECK(r_spher.radius == 1);
        CHECK(r_spher.phi == 0);
        CHECK(r_spher.theta == pi);
    }

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
    center_t r = {1,2,1};
    double theta = pi/4;
    double phi = 0;

    auto Y1 = eval_spherical_harmonics(q1,theta,phi);
    auto Y2 = eval_spherical_harmonics(q2,theta,phi);

    CHECK(abs(Y1 - std::complex<double>(-0.38627420,0)) < 10E-9) ;
    CHECK(abs(Y2 - std::complex<double>(0.38627420,0)) < 10E-9) ;
    {
        Spherical_Coordinates sc(r);
        auto Y3 = eval_spherical_harmonics(q1,sc.theta,sc.phi);
        CHECK(abs(Y3 - std::complex<double>(-0.128758, -0.257516)) < 10E-6);
    }

}

TEST_CASE( "Evaluate B Functions ", "[b_func_engine]" ) {
    ///TEST MAY NEED FIXING
    Quantum_Numbers q1 = {4,2,1};
    double alpha = 1.0;
    center_t r{1,2,1} ;
    {
        //check internal boost bessel function
        Spherical_Coordinates sc(r);
        CHECK(bm::cyl_bessel_k(q1.n - 1.0 / 2.0, alpha * sc.radius) - 0.48190552 < 10E-9);
    }
    B_function_Engine b_func_engine;
    auto B = b_func_engine.calculate(q1,alpha,r);
    CHECK(abs(B - std::complex<double>(-0.000148279, -0.000296558)) < 10E-9) ;

    {   //case where result is 0
        double alpha = 1;
        center_t r{1,1,0} ;
        B_function_Engine b_func_engine;
        auto B = b_func_engine.calculate(q1,alpha,r);
        CHECK(abs(B - std::complex<double>(0,0)) < 10E-9) ;
    }
    {   //case where r is 0
        double alpha = 1;
        center_t r{0,0,0} ;
        B_function_Engine b_func_engine;
        auto B = b_func_engine.calculate(q1,alpha,r);
        CHECK(abs(B - std::complex<double>(0,0)) < 10E-9) ;
    }/*
    {   //case where r is 0 , n=1
        Quantum_Numbers q2 = {1,2,1};
        double alpha = 1;
        center_t r{0,0,0};
        B_function_Engine b_func_engine;
        auto B = b_func_engine.calculate(q2,alpha,r);
        CHECK(abs(B - std::complex<double>(1,0)) < 10E-9) ;
    }*/

}
