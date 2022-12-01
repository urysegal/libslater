#pragma once
#include <vector>
#include "libslater.h"
#include "nested_summation.h"
#include "integrators.h"

namespace slater {

const std::string electron_repulsion_sdbar_name = "sdbar-electron-repulsion";

typedef int indexer_t; /// For this algorithm, the indices are integers (negative and positive )

class SDbar_Sum_State : public Summation_State<indexer_t, std::complex<double> > {

public:

    // The naming of the members here is done with the goal of one-to-one match with the paper in reference [4]
    // for easier reading and comparing.

    // Problem parameters
    std::array<unsigned int, 4> n_as_vec;
    std::array< int, 4> l_as_vec;
    std::array< int, 4> m_as_vec;

    unsigned int n1,n2,n3,n4;
    int l1,l2,l3,l4;
    int m1,m2,m3,m4;

    center_t A = {};
    center_t B = {};
    center_t C = {};
    center_t D = {};
    std::array<sto_exponent_t, 4> zeta;


    // Precalculated values

    center_t R21_vec;
    center_t R43_vec;
    double R21 ;
    double R34 ;


    // Indexing variables

    indexer_t l1_tag;
    indexer_t l2_tag;
    indexer_t l3_tag;
    indexer_t l4_tag;

    indexer_t m1_tag;
    indexer_t m2_tag;
    indexer_t m3_tag;
    indexer_t m4_tag;

    indexer_t l;
    indexer_t l12;
    indexer_t l_tag;
    indexer_t l34;

    auto get_m21() { return m2 - m2_tag - (m1 - m1_tag) ; }
    auto get_m43() { return m4 - m4_tag - (m3 - m3_tag) ; }



};

class Electron_Repulsion_SDbar : public STO_Integrator {

public:

    /// Build a 4-center 2-electron evaluator using SDbar transformation
    Electron_Repulsion_SDbar() ;

    /// This constructor is only used for building the factory
    explicit Electron_Repulsion_SDbar(const std::string &name) : STO_Integrator(name) {}

    /// Release any memory used by the integrator
    virtual ~Electron_Repulsion_SDbar() = default;

    [[nodiscard]] STO_Integrator *clone() const override;

    /// Initialize the integrator with a set of options
    /// \param params set of options for the integrator
    void init(const STO_Integration_Options &params) override;

    energy_unit_t integrate(
            const std::vector<STO_Basis_Function> &functions,
            const std::vector<center_t> &centers
    ) override;


private:

    SDbar_Sum_State state;
    int number_of_quadrature_points = 30; /// How many quadrature points we should calculate

    void setup_state(const std::vector<STO_Basis_Function> &functions);
    energy_unit_t evaluate();

};

}