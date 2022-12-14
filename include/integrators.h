#pragma once
#include "bfunctions.h"

namespace slater {
/// Implementation the integrals needed for HF or DFT calculation in STO basis set. A class implements
/// one of the required integrals. On template instantiation time, the user picks which one is implemented.
/// If there is a lot of common code betwween two integratal calculation, due to common technique (e.g. using
/// B-functions ), multiple inheritance can be used (or any other common-code technique ).

class STO_Integrator {

public:

    STO_Integrator(const std::string name_);

    virtual ~STO_Integrator() = default ;

    /// Initialize the integration engine with a set of (possibly empty) options
    /// \param options A set of option that may modify the behaviour of the algorithms in this class
    virtual void init(const STO_Integration_Options &options) = 0 ;

    /// Perform an integration calculation
    /// \param functions the functions to integrate over
    /// \param centers external points (e.g. nuclei position), if any
    /// \return the integral value
    virtual energy_unit_t integrate(
            const std::vector<STO_Basis_Function> &functions,
            const std::vector<center_t> &centers
    ) = 0;

    /// Return a new class from the factory
    /// \return Pointer to the STO integrator requested
    [[nodiscard]] virtual STO_Integrator *clone() const = 0;

    /// Given an integrator name, return a new instance of it
    /// \param name name of the integrator to create
    /// \return pointer to the integrator
    static STO_Integrator *create(const std::string &name);

    /// Create all the pair of B functions and their normalization coefficients from the two given sequences of B functions,
    /// each representing an STO. The result is kept in the "equivalence_series" member
    /// \param f1 First sequences of B functions
    /// \param f2 Second sequences of B functions
    void create_integration_pairs(const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2) ;


protected:
    const int number_of_centers = 0;
    const int number_of_electrons = 1;
    const int number_of_external_centers = 0;

    /// When integrating of each pair of functions in this vector and then summing up the values, you get the
    /// Overlap integral value
    std::vector<std::pair<
            std::pair<double,B_function_details>,
            std::pair<double, B_function_details> > > equivalence_series;


    static std::map<std::string, STO_Integrator *> all_integrators;

    /// Create an STO integrator of some kind
    /// \param number_of_centers_ how many different basis functions we have
    /// \param number_of_electrons_ how many electrons positions we integrate over
    /// \param number_of_external_centers_ how many external centers are involved.
    STO_Integrator( int number_of_centers_, int number_of_electrons_ = 1, int number_of_external_centers_ = 0 ) :

            number_of_centers(number_of_centers_), number_of_electrons(number_of_electrons_),
            number_of_external_centers(number_of_external_centers_)
    {
    }

    /// Do an integration of two STO using their B functions representation.
    /// \param f1 first representation of STO by B functions
    /// \param f2 second representation of STO by B functions
    /// \return  integral value
    std::complex<double> do_integrate( const B_functions_representation_of_STO &f1, const B_functions_representation_of_STO &f2 ) ;

    /// Do an integration using the two B functions given. The specific integral to calculate is defined on the subclass.
    /// \param f1 first B function
    /// \param f2 second B functions
    /// \return  integral value
    [[nodiscard]] virtual std::complex<double> integrate_using_b_functions(const B_function_details &f1, const B_function_details &f2)  = 0;


};

}