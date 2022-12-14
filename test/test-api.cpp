#include <libslater.h>
#include <math.h>
#include <gaunt.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

using namespace slater;

class testRunListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void testRunStarting(Catch::TestRunInfo const&) override {
        slater::libslater_global_init();
        srand(time(NULL));
    }

    void testRunEnded( Catch::TestRunStats const& testRunStats ) override {
        slater::libslater_global_cleanup();
    }

};

CATCH_REGISTER_LISTENER(testRunListener)

TEST_CASE( "One overlap integral", "[overlap]" ) {

    Quantum_Numbers quantum_numbers = {2,1,0};

    STO_Basis_Function_Info oxygen_s( 0.252, quantum_numbers);
    STO_Basis_Function_Info hydrogen_s( 0.952, quantum_numbers);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});


    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    if (engine) {

        STO_Integration_Options parameters;
        parameters.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(parameters);

        energy_unit_t result = engine->overlap({oxygen_1_s, hydrogen_1_s});
        CHECK(result.imag() == 0 );
        CHECK(abs(result.real() - 0.3400571077) < 10e-7 );

        delete engine;

    }

}

TEST_CASE( "One nuclear attraction integral", "[nuclear]" ) {

    Quantum_Numbers quantum_numbers = {2,1,0};

    STO_Basis_Function_Info oxygen_s( 0.252, quantum_numbers);
    STO_Basis_Function_Info hydrogen_s( 0.952, quantum_numbers);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});


    STO_Integration_Engine engine_factory;
    std::map<slater::integration_types, std::string> engines;
    auto engine = engine_factory.create(engines);
    if (engine) {

        STO_Integration_Options parameters;
        parameters.set(Number_of_quadrature_points_Parameter_Name,30);
        engine->init(parameters);

        auto result = engine->nuclear_attraction({oxygen_1_s, hydrogen_1_s}, {0,0.5,-1});
        CHECK(result.imag() == 0 );
        //CHECK(abs(result.real() + 93.2411967435 ) < 0.000001 );
        //CHECK(abs(result.real() == 0 );
        delete engine;
    }
}


TEST_CASE( "nonexistent engine", "[api]" )
{
    std::map<slater::integration_types, std::string> engines;
    engines.emplace(slater::integration_types::OVERLAP, "blah blah");
    auto engine = STO_Integration_Engine().create(engines);
    CHECK(engine == nullptr);
}

TEST_CASE("Options behavior", "[api]")
{
    STO_Integration_Options options;

    bool no_value = false;
    CHECK(options.get("no such option as this", no_value) == false);

    {
        options.set(Use_Normalized_B_Functions_Parameter_Name, true);
        bool check_value = false;
        options.get(Use_Normalized_B_Functions_Parameter_Name, check_value);
        CHECK(check_value == true);
    }

    {
        options.set(Number_of_quadrature_points_Parameter_Name, int(180));
        int check_value = false;
        options.get(Number_of_quadrature_points_Parameter_Name, check_value);
        CHECK(check_value == 180);
    }

    {
        options.set(Number_of_quadrature_points_Parameter_Name, (double)0.45);
        double check_value = false;
        options.get(Number_of_quadrature_points_Parameter_Name, check_value);
        CHECK(check_value == 0.45);
    }

}

TEST_CASE("Getters and Setters", "[api]")
{
    Quantum_Numbers quantum_numbers = {2,1,0};
    STO_Basis_Function_Info oxygen_s(0.252, quantum_numbers);
    oxygen_s.set_coefficient(1.0);
    oxygen_s.set_exponent(2.0);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});
    CHECK(oxygen_1_s.get_normalization_coefficient() == 1);
    CHECK(oxygen_1_s.get_exponent() == 2);

    Quantum_Numbers nqn = { 3,2,1, spin_quantum_number_t::UP};
    oxygen_s.set_quantum_numbers(nqn);
    auto qn = oxygen_s.get_quantum_numbers();
    CHECK(qn.n == 3);
    CHECK(qn.l == 2);
    CHECK(qn.m == 1);
    CHECK(qn.ms == spin_quantum_number_t::UP);

}

TEST_CASE("Gaunt Coefficients", "[api]")
{
    auto g = slater::Gaunt_Coefficient_Engine::get()->calculate({1,1,0,0,1,-1});
    CHECK(g-0.282095 < 0.00001);
    g = slater::Gaunt_Coefficient_Engine::get()->calculate({6,0,6,6,6,-6});
    CHECK(g- -0.0629787762403749 < 0.00000001);
    g = slater::Gaunt_Coefficient_Engine::get()->calculate({5,4,3,-1,2,-3});
    CHECK(g < 0.00000001);

}

TEST_CASE( "evaliate STO function", "[STO]" ) {

    Quantum_Numbers quantum_numbers = {2, 1, 0};

    STO_Basis_Function_Info oxygen_s(0.252, quantum_numbers);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});

    CHECK(abs(oxygen_1_s.evaluate({1,1,1}).real() - 0.0116243037) < 10e-9 );
    CHECK(abs(oxygen_1_s.evaluate_conjugate({1,1,1}).real() - 0.0116243037) < 10e-9);
    CHECK(oxygen_1_s.evaluate({1,1,1}).imag() == 0 );
    CHECK(oxygen_1_s.evaluate_conjugate({1,1,1}).imag() == 0 );


}