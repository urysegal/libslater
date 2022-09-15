#include <libslater.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/interfaces/catch_interfaces_reporter.hpp>

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

TEST_CASE( "bad call to API", "[overlap]" ) {
}

