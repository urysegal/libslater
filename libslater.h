#pragma once
#include <vector>
#include <array>
#include <string>



namespace slater {

/// Parameter to specify if we use normalized B functions
static constexpr const char *Use_Normalized_B_Functions_Parameter_Name = "use_normalized_b_functions" ;
/// Parameter to specify number of quadrature points
static constexpr const char *Number_of_quadrature_points_Parameter_Name = "number_of_quadrature_points" ;


/// Initialize the library
void libslater_global_init();

/// Clean up the library
void libslater_global_cleanup();

/// Some quantity of energy, which is returned from the integration routines
typedef double energy_unit_t;

/// An STO exponent
typedef double sto_exponent_t;

/// A coefficient or Normalization coefficient
typedef double sto_coefficient_t;

typedef unsigned int angular_quantum_number_t;
typedef unsigned int principal_quantum_number_t;
typedef int magnetic_quantum_number_t;

/// A coordinate in 1D space
typedef double spatial_coordinate_t;

/// A coordinate in 3D space
typedef std::array<spatial_coordinate_t, 3> center_t;

enum class spin_quantum_number_t  { UNDEFINED, UP, DOWN } ;


struct Quantum_Numbers {
    principal_quantum_number_t n = 0; /// principal quantum number, also known as 'n'
    angular_quantum_number_t l = 0; /// Angular moment number, also known as 'l'
    magnetic_quantum_number_t m = 0; /// Magnetic/Orientation quantum number, also known as 'm' or 'ml'
    spin_quantum_number_t ms = spin_quantum_number_t::UNDEFINED; /// Spin Quantum number, a.k.a. 'ms'
};

/// One STO style basis function details. STO have radial part and angular part.
/// the radial part is the form f(r)=N*r^(n-1)*exp(-ar). N and a are the first two parameters of the constructor.
/// the angular part is a spherical harmonic which is given by the third parameter.

class STO_Basis_Function_Info {

    sto_coefficient_t coefficient; /// Coefficient of the radial part. also known as N
    sto_exponent_t exponent; /// exponent of the radial part. also known as alpha
    Quantum_Numbers quantum_numbers;     /// Quantum numbers for this basis function

public:

    /// Get the set of quantum numbers for this basis function
    /// \return set of quantum numbers for this basis function
    const Quantum_Numbers &get_quantum_numbers() const;

    /// Set the quantum numbers for this basis function
    /// \param quantum_numbers set of quantum numbers to use
    void set_quantum_numbers(const Quantum_Numbers &quantum_numbers);

    /// Get the exponent used by this basis function
    /// \return the exponent used by this basis function
    sto_exponent_t get_exponent() const ;

    /// Set the exponent used by this basis function
    /// \param e exponent to use
    void set_exponent(sto_exponent_t e) ;

    /// Get the coefficient used by this basis function
    /// \return the coefficient used by this basis function
    sto_coefficient_t get_coefficient() const ;

    /// Set the coefficient used by this basis function
    /// \param c coefficient to use
    void set_coefficient(sto_coefficient_t c) ;


public:
    /// Construct an STO style basis function detail object.
    /// \param coefficient_  Coefficient of the radial part.
    /// \param exponent_ exponent of the radial part.
    /// \param quantum_numbers_ Set of quantum numbers for this function
    STO_Basis_Function_Info(sto_coefficient_t coefficient_, sto_exponent_t exponent_, const Quantum_Numbers &quantum_numbers);

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

    /// Get the set of quantum numbers for this basis function
    /// \return set of quantum numbers
    const Quantum_Numbers &get_quantum_numbers() const ;

    /// Get the exponent used by this basis function
    /// \return the exponent used by this basis function
    sto_exponent_t get_exponent() const ;

    /// Get the coefficient used by this basis function
    /// \return the coefficient used by this basis function
    sto_coefficient_t get_coefficient() const ;

    /// Get the spatial center of this basis function
    /// \return the spatial center of this basis function
    center_t get_center() const ;


};

class STO_Integration_Options_Impl;

/// Class to maintain the various parameters and options that can be used to twick integration implementations
class STO_Integration_Options
{

    STO_Integration_Options_Impl *implementation = nullptr; /// Implementation-dependant storage for option values

public:

    STO_Integration_Options();

    virtual ~STO_Integration_Options();

    /// set ( or override previous setting )  a parameter
    /// \param name parameter name to set
    /// \param value  value to set
    template <class T> void set(const std::string &name, const T& value);

    /// Get a boolean parameter by name. If the parameter is not set, "value" will not be touched, so that
    /// default values can be kept as it.
    /// \param name parameter name to get
    /// \param value parameter value returned in this reference
    /// \return true if this parameter was set at all. If false, there is no value given to "name"
    template <class T> bool get(const std::string &name, T &value) const;
};


/// Implementation of all the integrals needed for HF or DFT calculation in STO basis set.

class STO_Integrator {

public:
    STO_Integrator() = default;
    virtual ~STO_Integrator() = default ;

    /// Initialize the integration engine with a set of (possibly empty) options
    /// \param options A set of option that may modify the behaviour of the algorithms in this class
    virtual void init(const STO_Integration_Options &options) = 0 ;

    /// Calculate the Overlap Integral <f|g> over the given two STO basis functions
    /// \param functions The two function whose overlap is to be calculated
    /// \return The value of the overlap integral
    virtual energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &functions) = 0;
};


/// Integrator creation factory that creates the desired type of STO integrator. Don't forget to free the pointer returned by create()
class STO_Integration_Engine {

public:
    /// Construct an integration engines factory
    STO_Integration_Engine() ;

    /// Create an STO integration engine. Use "default" as engine type to get the default implementation
    /// or your preferred implementation if you have another.
    /// \param engine_type name of engine to create. Use "default" for the library default one (currently Homeier B-functions)
    /// \return a pointer to the engine, which must be released by called to free resources. nullptr if there was no implementation matching the string given
    STO_Integrator *create(const std::string &engine_type);
};

}