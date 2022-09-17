#pragma once
#include <vector>
#include <array>
#include <string>

namespace slater {

/// Initialize the library
void libslater_global_init();

/// Clean up the library
void libslater_global_cleanup();

typedef double sto_exponent_t;
typedef double sto_coefficient_t;
typedef unsigned int moment_quantum_number_t;

typedef double spatial_coordinate_t;
typedef std::array<spatial_coordinate_t, 3> center_t;

typedef double integral_value;


/// One STO style basis function details. STO have radial part and angular part.
/// the radial part is the form f(r)=N*r^(n-1)*exp(-ar). N and a are the first two parameters of the constructor.
/// the angular part is a spherical harmonic which is given by the third parameter.
/// n is limited to be equal to l in this library, as in GTO-based basis sets.

class STO_Basis_Function_Info {

    sto_coefficient_t coefficient; /// Coefficient of the radial part.
    sto_exponent_t exponent; /// exponent of the radial part.
    moment_quantum_number_t moment_quantum_number; /// Angular moment number, also known as 'l'

public:
    /// Contruct an STO style basis function detail object.
    /// \param coefficient_  Coefficient of the radial part.
    /// \param exponent_ exponent of the radial part.
    /// \param moment_quantum_number_ Angular moment number
    STO_Basis_Function_Info(sto_coefficient_t coefficient_, sto_exponent_t exponent_, moment_quantum_number_t moment_quantum_number_);

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

class STO_Integration_Parameters
{
public:
    void add(const std::string &name, bool value);
};


class STO_Integrator {

public:
    STO_Integrator();
    virtual ~STO_Integrator() = default;

    virtual void init(const STO_Integration_Parameters &params) = 0 ;
    virtual integral_value overlap(const std::array<STO_Basis_Function, 4> &) = 0;
};


class STO_Integration_Engine {

public:
    STO_Integration_Engine();
    STO_Integrator *create(const std::string &engine_type);

};




}

