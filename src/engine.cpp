#include "libslater.h"
#include "homeier.h"
#include "logger.h"

namespace slater {

const std::string default_engine = "B-functions-homeier";


STO_Integration_Engine::STO_Integration_Engine()
{
    logger()->debug("Engine factory created");
}

STO_Integrator *STO_Integration_Engine::create(const std::string &engine_type)
{
    STO_Integrator *res = nullptr;
    std::string type_to_use = engine_type;
    if ( engine_type == "default" ) {
        type_to_use = default_engine;
    }
    if ( type_to_use == "B-functions-homeier" ) {
        res = new Homeier_Integrator();
    }

    if ( res ) {
        logger()->info("{}-type engine was created", type_to_use);
    } else {
        logger()->info("Could not find {}-type engine", type_to_use);
    }

    return res;
}

}
