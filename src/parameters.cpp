#include <map>
#include <any>
#include "libslater.h"



namespace slater {


class STO_Integration_Options_Impl
{

    std::map<std::string , std::any> all_options;
public:

    STO_Integration_Options_Impl() = default ;
    virtual ~STO_Integration_Options_Impl() = default;


    template<class T> void set(const std::string &name, const T& value)
    {
        all_options[name] = value;
    }


    template<class T> bool get(const std::string &name, T &value) const
    {
        bool res = true;
        auto it = all_options.find(name);
        if ( it == all_options.end() ) {
            res = false;
        } else {
            value = any_cast<T>(it->second);
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

template <class T>
void STO_Integration_Options::set(const std::string &name, const T& value)
{
    implementation->set(name, value);
}

template <class T>
    bool STO_Integration_Options::get(const std::string &name, T &value) const
{
    return implementation->get(name, value);
}

}
