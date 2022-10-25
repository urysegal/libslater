#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"

#include <iomanip>

using namespace slater;

int
main()
{


    Quantum_Numbers quantum_numbers = {2,1,0};

    STO_Basis_Function_Info oxygen_s( 1, quantum_numbers);
    STO_Basis_Function_Info oxygen_p(2, quantum_numbers);

    STO_Basis_Function_Info hydrogen_s( 1, quantum_numbers);


    STO_Basis_Function oxygen_1_s(oxygen_s, {2, 0, 0});
    STO_Basis_Function oxygen_2_p(oxygen_p, {0, 0, 1});
    STO_Basis_Function hydrogen_1_s(hydrogen_s, {0.70710678, 0, 0.56568542});
    STO_Basis_Function hydrogen_2_s(hydrogen_s, {-0.70710678, 0, 0.56568542});

    std::vector<STO_Basis_Function> basis_set = {oxygen_1_s, oxygen_2_p, hydrogen_1_s, hydrogen_2_s};

    std::map<slater::integration_types, std::string> engines;
    auto engine = STO_Integration_Engine().create(engines);
    if (engine) {

        STO_Integration_Options options;
        options.set(Use_Normalized_B_Functions_Parameter_Name, true);

        engine->init(options);
       auto result = engine->overlap({basis_set[0], basis_set[1]});


        std::cout << std::fixed;
        std::cout << std::setprecision(16);

       std::cout << result << std::endl;
        //const double order = 5;
       // const double z[3] ={ 2.2532185649430092, 7.3753971484553746 , 22.190479333842703 };

      //  for ( int i = 0 ; i < 3 ; ++i ) {
        //    do_compute_reduced_bessel_function_half(
          //          order,
      //              z[i]
        //    );


            //printf("Bessel %15.15f\n", r);
        //}
      //  calculate_gauss_point(basis_set[0], basis_set[1], 1.5554323515926116E-004);

        //result = engine->nuclear_attraction({basis_set[2], basis_set[0]}, {1,0,0.5});
        //std::cout << result << std::endl;

        delete engine;

    }

	return 0;
}
