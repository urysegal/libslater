#include "libslater.h"
#include "homeier.h"
#include "analytical-3c.h"
#include "logger.h"

namespace slater {

std::map<std::string, STO_Integrator *> STO_Integrator::all_integrators;

const std::map<integration_types, std::string > default_engines =
        {
                {integration_types::OVERLAP, overlap_homeier_imp_name},
                {integration_types::KINETIC, kinetic_homeier_imp_name},
                {integration_types::NUCLEAR_ATTRACTION, analytical_3c_name }
        };

STO_Integrator *STO_Integrator::create(const std::string &name)
{
    STO_Integrator *res = nullptr;

    auto it = all_integrators.find(name);
    if ( it != all_integrators.end() ) {
        res = it->second->clone();
    }
    return res;
}


void STO_Integrations::add_engine(slater::integration_types type, slater::STO_Integrator *integrator)
{
    this->integrators.emplace(type, integrator);
}

STO_Integration_Engine::STO_Integration_Engine()
{
    logger()->debug("Engine factory created");
}


STO_Integrator::STO_Integrator(const std::string name_)
{
    all_integrators.emplace(name_, this);
}

void STO_Integrator::create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2)
{
    equivalence_series.clear();
    for ( auto i : f1) {
        for (auto j: f2) {
            equivalence_series.emplace_back(i, j);
        }
    }
}

complex STO_Integrator::do_integrate( const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2 )
{
    create_integration_pairs(f1, f2);

    std::vector<energy_unit_t> partial_results;
    for (auto const &p: equivalence_series) {
        // int(B_1*B_2)
        energy_unit_t partial_result = integrate_using_b_functions(p.first.second,p.second.second);

        // Bcoeff1*Bcoeff2 * int(B_1*B_2)
        partial_result *= p.first.first * p.second.first;
        partial_results.emplace_back(partial_result);
    }

    // sum(sum(Bcoeff1*Bcoeff2 * int(B_1 B_2)))
    energy_unit_t result = 0;
    for (auto &pr: partial_results)
        result += pr;

    result *= f1.get_rescaling_coefficient() * f2.get_rescaling_coefficient();

    return result;

}

STO_Integrations *STO_Integration_Engine::create(std::map<integration_types, std::string > &engines)
{
    auto res = new STO_Integrations;

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


energy_unit_t STO_Integrations::two_functions_integral(const std::array<STO_Basis_Function, 2> &functions,
                                                       const center_t &nuclei,
                                                       integration_types which_type)
{

    assert(functions[0].get_quantum_numbers().l + functions[1].get_quantum_numbers().l
           <= Gaunt_Coefficient_Engine::get_maximal_gaunt_l());

    functions[0].get_quantum_numbers().validate();
    functions[1].get_quantum_numbers().validate();


    auto it = this->integrators.find(which_type);
    if ( it != this->integrators.end() ) {
        return it->second->integrate({functions[0], functions[1]}, {nuclei});
    } else {
        throw std::runtime_error("Cannot find integral implementation"); // LCOV_EXCL_LINE
    }
}


energy_unit_t STO_Integrations::overlap(const std::array<STO_Basis_Function, 2> &functions)
{
    return two_functions_integral(functions, {},integration_types::OVERLAP);
}


energy_unit_t STO_Integrations::kinetic(const std::array<STO_Basis_Function, 2> &functions)
{
    return two_functions_integral(functions, {},integration_types::KINETIC);
}


energy_unit_t STO_Integrations::nuclear_attraction(const std::array<STO_Basis_Function, 2> &functions,
                                                   const center_t &nuclei)
{
    return two_functions_integral(functions, nuclei,integration_types::NUCLEAR_ATTRACTION);
}


void STO_Integrations::init(const slater::STO_Integration_Options &options) {
    for (auto it: this->integrators) {
        it.second->init(options);
    }
}



}
