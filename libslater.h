#pragma once
#include <vector>
#include <array>
#include <string>

namespace slater {

/// Parameter to specify if we use normalized B functions
static constexpr const char *Use_Normalized_B_Functions_Parameter_Name = "use_normalized_b_functions" ;


/// Initialize the library
void libslater_global_init();

/// Clean up the library
void libslater_global_cleanup();

typedef double sto_exponent_t;
typedef double sto_coefficient_t;

typedef unsigned int moment_quantum_number_t;
typedef int principal_quantum_number_t;
typedef int orientation_quantum_number_t;

typedef double spatial_coordinate_t;
typedef std::array<spatial_coordinate_t, 3> center_t;

typedef double integral_value;


/// One STO style basis function details. STO have radial part and angular part.
/// the radial part is the form f(r)=N*r^(n-1)*exp(-ar). N and a are the first two parameters of the constructor.
/// the angular part is a spherical harmonic which is given by the third parameter.

class STO_Basis_Function_Info {

    sto_coefficient_t coefficient; /// Coefficient of the radial part. also known as N
    sto_exponent_t exponent; /// exponent of the radial part.
    moment_quantum_number_t moment_quantum_number; /// Angular moment number, also known as 'l'
    principal_quantum_number_t principal_quantum_number; /// principal quantum number, also known as 'n'
    orientation_quantum_number_t orientation; /// also known as 'm'
public:
    /// Contruct an STO style basis function detail object.
    /// \param coefficient_  Coefficient of the radial part.
    /// \param exponent_ exponent of the radial part.
    /// \param moment_quantum_number_ Angular moment number, also known as 'l'
    /// \param principal_quantum_number_ principal quantum number, also known as 'n'
    /// \param orientation orientation of the function, also known as 'm'
    STO_Basis_Function_Info(sto_coefficient_t coefficient_, sto_exponent_t exponent_, moment_quantum_number_t moment_quantum_number_,
                            principal_quantum_number_t principal_quantum_number_, orientation_quantum_number_t orientation_);

};


/// One STO basis function, located at a given center.
class STO_Basis_Function  {

    STO_Basis_Function_Info function_info; /// Basis function parameters
    center_t center; /// the center of the function

public:

    /// Construct a basis function located at a specific coordinates
    /// \param function_info Basis function information
    /// \param location Cartesian center of the function
    STO_Basis_Function(STO_Basis_Function_Info function_info_, center_t location_);

};

class STO_Integration_Options_Impl;

/// Class to maintain the various parameters and options that can be used to twick integration implementations
class STO_Integration_Options
{

    STO_Integration_Options_Impl *implementation = nullptr;

public:

    STO_Integration_Options();

    virtual ~STO_Integration_Options();

    /// set ( or override previous setting ) of a boolean parameter
    /// \param name parameter name to set
    /// \param value paramter boolean value to set
    virtual void set(const std::string &name, bool value);

    /// Get a boolean parameter by name. If the parameter is not set, "value" will not be touched, so that
    /// default values can be kept as it.
    /// \param name parameter name to get
    /// \param value parameter value returned in this reference
    /// \return true if this parameter was set at all. If false, there is no value given to "name"
    virtual bool get(const std::string &name, bool &value) const;
};


/// Implementation of all the integrals needed for HF or DFT calculation in STO basis set.

class STO_Integrator {

public:
    STO_Integrator();
    virtual ~STO_Integrator() ;

    /// Initialize the integration engine with a set of (possibly empty) options
    /// \param options A set of option that may modify the behaviour of the algorithms in this class
    virtual void init(const STO_Integration_Options &options) = 0 ;

    /// Calculate the Overlap Integral <f|g> over the given two STO basis functions
    /// \param functions The two funciton whose overlap is to be calculated
    /// \return The value of the overlap integral
    virtual integral_value overlap(const std::array<STO_Basis_Function, 2> &functions) = 0;
};


/// Integrator creation factory that creates the desired type of STO integrator. Don't forget to free the pointer retunred by create()
class STO_Integration_Engine {

public:
    STO_Integration_Engine();

    /// Create an STO integration engine. Use "default" as engine type to get the default implementation
    /// or your preferred implementation if you have another.
    /// \param engine_type name of engine to create. Use "default" for the library default one (currently Homeier B-functions)
    /// \return a pointer to the engine, which must be released by called to free resources. nullptr if there was no implementation matching the string given
    STO_Integrator *create(const std::string &engine_type);
};

}

