#pragma once
#include <vector>
#include <array>
#include <string>
#include <complex>
#include <map>



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
typedef std::complex<double> energy_unit_t;

/// An STO alpha
typedef double sto_exponent_t;

/// A normalization_coefficient or Normalization normalization_coefficient
typedef double sto_coefficient_t;

typedef  int angular_quantum_number_t;
typedef  int principal_quantum_number_t;
typedef int magnetic_quantum_number_t;

/// A coordinate in 1D space
typedef double spatial_coordinate_t;

/// A coordinate in 3D space
typedef std::array<spatial_coordinate_t, 3> center_t;

enum class spin_quantum_number_t  { UNDEFINED, UP, DOWN } ;

enum class integration_types : int {
    OVERLAP,
    KINETIC,
    NUCLEAR_ATTRACTION
};


struct Quantum_Numbers {
    principal_quantum_number_t n = 0; /// principal quantum number, also known as 'n'
    angular_quantum_number_t l = 0; /// Angular moment number, also known as 'l'
    magnetic_quantum_number_t m = 0; /// Magnetic/Orientation quantum number, also known as 'm' or 'ml'
    spin_quantum_number_t ms = spin_quantum_number_t::UNDEFINED; /// Spin Quantum number, a.k.a. 'ms'
    void validate() const;
};

/// One STO style basis function details. STO have radial part and angular part.
/// the radial part is the form f(r)=N*r^(n-1)*exp(-ar). N and a are the first two parameters of the constructor.
/// the angular part is a spherical harmonic which is given by the third parameter.

class STO_Basis_Function_Info {

    sto_coefficient_t normalization_coefficient; /// Normalization Coefficient N(n,alpha) of the radial part. also known as N
    sto_exponent_t exponent; /// alpha of the radial part. also known as alpha
    Quantum_Numbers quantum_numbers;     /// Quantum numbers for this basis function

public:

    /// Get the set of quantum numbers for this basis function
    /// \return set of quantum numbers for this basis function
    const Quantum_Numbers &get_quantum_numbers() const;

    /// Set the quantum numbers for this basis function
    /// \param quantum_numbers set of quantum numbers to use
    void set_quantum_numbers(const Quantum_Numbers &quantum_numbers);

    /// Get the alpha used by this basis function
    /// \return the alpha used by this basis function
    sto_exponent_t get_exponent() const ;

    /// Set the alpha used by this basis function
    /// \param e alpha to use
    void set_exponent(sto_exponent_t e) ;

    /// Get the normalization_coefficient used by this basis function
    /// \return the normalization_coefficient used by this basis function
    sto_coefficient_t get_coefficient() const ;

    /// Set the normalization_coefficient used by this basis function
    /// \param c normalization_coefficient to use
    void set_coefficient(sto_coefficient_t c) ;


public:
    /// Construct an STO style basis function detail object.
    /// \param exponent_ alpha of the radial part.
    /// \param quantum_numbers Set of quantum numbers for this function
    STO_Basis_Function_Info( sto_exponent_t exponent_, const Quantum_Numbers &quantum_numbers);

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

    /// Get the alpha used by this basis function
    /// \return the alpha used by this basis function
    sto_exponent_t get_exponent() const ;

    /// Get the normalization_coefficient used by this basis function
    /// \return the normalization_coefficient used by this basis function
    sto_coefficient_t get_coefficient() const ;

    /// Get the spatial center of this basis function
    /// \return the spatial center of this basis function
    center_t get_center() const ;

    /// Return the value of the function at the given point
    /// \param point point to evaliate at
    /// \return complex value of function at the given point
    std::complex<double> evaluate(const center_t &point) const;

    /// Return the value of the conjugate function at the given point
    /// \param point point to evaliate at
    /// \return complex value of function at the given point
    std::complex<double> evaluate_conjugate(const center_t &point) const;


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

class STO_Integrator;

/// Provide implementations of all integral calculation needed for Grid-Free DFT and HF.
class STO_Integrations {

public:
    STO_Integrations() = default;
    ~STO_Integrations() ;

    /// Initialize the integration engine with a set of (possibly empty) options
    /// \param options A set of option that may modify the behaviour of the algorithms in this class
    void init(const STO_Integration_Options &options) ;

    void add_engine(integration_types type, STO_Integrator *integrator);


    /// Calculate the Overlap Integral <f|g> over the given two STO basis functions
    /// \param functions The two function whose overlap is to be calculated
    /// \return The value of the overlap integral
    energy_unit_t overlap(const std::array<STO_Basis_Function, 2> &functions) ;

    /// Calculate the Kinetic Energy Integral <f|del^2|g> over the given two STO basis functions
    /// \param functions The two function whose overlap is to be calculated
    /// \return The value of the kinetic energy integral
    energy_unit_t kinetic(const std::array<STO_Basis_Function, 2> &functions) ;

    /// Calculate the nuclear attraction integral <f|1/R|g> over the given two STO basis functions and the one nuclei
    /// \param functions The two function whose overlap is to be calculated
    /// \param nuclei the position of the nuclei
    /// \return The value of the nuclear attraction integral
    energy_unit_t nuclear_attraction(const std::array<STO_Basis_Function, 2> &functions, const center_t &nuclei) ;

private:

    std::map<integration_types, STO_Integrator * > integrators;

};




/// Integrator creation factory that creates the desired type of STO integrator. Don't forget to free the pointer returned by create()
class STO_Integration_Engine {

public:


    /// Construct an integration engines factory
    STO_Integration_Engine() ;

    /// Create an STO integration engine. Use an empty map if you want to use the default engines, or pick specific
    /// implementation for specific integrals you'd like to override.
    /// \param engines a map from an engine type to the desired implementaion. empty map will use the default engines.
    /// \return a pointer to the engine, which must be released by called to free resources. nullptr if there was no
    /// implementation matching one or more of the strings given.
    STO_Integrations *create(std::map<integration_types, std::string > &engines);
};

}