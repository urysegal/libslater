#include "libslater.h"
#include "homeier.h"

namespace slater {

const std::string default_engine = "B-functions-homeier";


STO_Integration_Engine::STO_Integration_Engine()
{
}


STO_Integrator::STO_Integrator()
{
}

STO_Integrator *STO_Integration_Engine::create(const std::string &engine_type)
{
    std::string type_to_use =engine_type;
    if ( engine_type == "default" ) {
        type_to_use = default_engine;
    }
    if ( engine_type == "B-functions-homeier" ) {
        return new Homeier_Integrator();
    }
    return nullptr;
}

}
