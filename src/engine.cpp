#include "libslater.h"
#include "homeier.h"
#include "logger.h"

namespace slater {



const std::map<integration_types, std::string >  default_engines =
        {
                {integration_types::OVERLAP, overlap_homeier_imp_name},
        };


STO_Integration_Engine::STO_Integration_Engine()
{
    logger()->debug("Engine factory created");
}


STO_Integrator::STO_Integrator(const std::string name_)
{
    this->all_integrators.emplace(name_, this);
}

STO_Integrations *STO_Integration_Engine::create(std::map<integration_types, std::string > &engines)
{
    STO_Integrations *res = new STO_Integrations;

    for ( auto &it : default_engines )
    {
        auto engine_type =it.first;
        auto engine_imp_name = it.second;

        auto alternative = engines.find(engine_type);
        if ( alternative != engines.end() ) {
            engine_imp_name = alternative->second;
        }

        auto engine_imp = STO_Integrator::create(engine_imp_name);
        if ( engine_imp )
        {
            logger()->info("{}-type integrator was created", engine_imp_name);
            res->add_engine(engine_type, engine_imp );
        } else {
            logger()->info("Could not find {}-type engine", engine_imp_name);
            delete res;
            res = nullptr;
            break;
        }

    }

    return res;
}


STO_Integrations::~STO_Integrations()
{
    for ( auto it: this->integrators ) {
        delete it.second;
    }
}

energy_unit_t STO_Integrations::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    auto it = this->integrators.find(integration_types::OVERLAP);
    if ( it != this->integrators.end() ) {
        return it->second->integrate({functions[0], functions[1]}, {});
    } else {
        throw std::runtime_error("Cannot find overlap integral implementation");
    }
}

void STO_Integrations::init(const slater::STO_Integration_Options &options) {
    for (auto it: this->integrators) {
        it.second->init(options);
    }
}

}
