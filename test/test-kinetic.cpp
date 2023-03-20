#include <libslater.h>
#include <math.h>
#include <gaunt.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>
#include "homeier.h"
using namespace slater;

double epsilon = 1e-8;

TEST_CASE( "One Kinetic B function integral", "[kinetic]" ) {

    Homeier_Integrator integrator(integration_types::KINETIC);
    const Quantum_Numbers q1 = {1,0,0};
    const Quantum_Numbers q2 = {1,0,0};

    center_t c1 = {2, 0, 1};
    center_t c2 = {0, 0, 1};

    B_function_details f1(1.0, q1, c1);
    B_function_details f2(1.0, q2, c2);

    auto res = 0;
    CHECK(res == 0.0);
}

