#include <map>
#include "libslater.h"



namespace slater {


class STO_Integration_Options_Impl
{

    std::map<std::string, bool> bool_options;

public:

    STO_Integration_Options_Impl() = default ;
    virtual ~STO_Integration_Options_Impl() = default;

    void set(const std::string &name, bool value)
    {
        bool_options[name] = value;
    }

    bool get(const std::string &name, bool &value) const
    {
        bool res = true;
        auto it = bool_options.find(name);
        if ( it == bool_options.end() ) {
            res = false;
        } else {
            value = it->second;
        }
        return res;
    }
};



STO_Integration_Options::STO_Integration_Options()
{
    implementation = new STO_Integration_Options_Impl();
}


STO_Integration_Options::~STO_Integration_Options()
{
    delete implementation;
}

void STO_Integration_Options::set(const std::string &name, bool value)
{
    implementation->set(name, value);
}

bool STO_Integration_Options::get(const std::string &name, bool &value) const
{
    return implementation->get(name, value);
}

}
