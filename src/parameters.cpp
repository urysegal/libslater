#include <map>
#include "libslater.h"



namespace slater {


class STO_Integration_Options_Impl : public STO_Integration_Options
{

    std::map<std::string, bool> bool_options;

public:

    STO_Integration_Options_Impl() = default ;
    virtual ~STO_Integration_Options_Impl() = default;

    virtual void set(const std::string &name, bool value) override
    {
        bool_options[name] = value;
    }

    virtual bool get(const std::string &name, bool &value) const override
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

}
