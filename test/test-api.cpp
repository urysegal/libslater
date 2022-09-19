#include <libslater.h>

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
    STO_Basis_Function_Info oxygen_s(43.5, 0.252, 0);
    STO_Basis_Function_Info hydrogen_s(3.15, 0.952, 0);

    STO_Basis_Function oxygen_1_s(oxygen_s, {0, 0, -0.14142136});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});


    auto engine = STO_Integration_Engine().create("default");
    if (engine) {

        STO_Integration_Options parameters;
        parameters.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(parameters);
        integral_value result = engine->overlap({ oxygen_1_s, hydrogen_1_s });
        delete engine;
        CHECK(result == 0);
    }

}

TEST_CASE( "nonexistent engine", "[api]" )
{
    auto engine = STO_Integration_Engine().create("blah blah");
    CHECK(engine == nullptr);
}

TEST_CASE("Options behavior", "[api]")
{
    STO_Integration_Options options;
    options.set(Use_Normalized_B_Functions_Parameter_Name, true);
    bool check_value = false;
    options.get(Use_Normalized_B_Functions_Parameter_Name, check_value);
    CHECK(check_value == true );
    CHECK(options.get("no such option as this", check_value) == false );
}