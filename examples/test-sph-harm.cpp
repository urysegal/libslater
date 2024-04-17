#include <iostream>
#include <Eigen/Dense>
#include "libslater.h"
#include <iomanip>
#include "libslater.h"
#include "coordinates.h"
#include "slater-utils.h"
using namespace slater;


int
main()
{
    principal_quantum_number_t n = 2;
    angular_quantum_number_t l = 2;
    std::complex<double> Yp;
    std::complex<double> Yn;
    for (int m=-l; m<=l; m++) {
        Quantum_Numbers q1 = {n, l, m};
        Quantum_Numbers q2 = {n, l, -1*m};
        center_t cart = {1,1, 0};
        Spherical_Coordinates sph(cart);

        Yp = eval_spherical_harmonics(q1, sph);
        Yn = eval_spherical_harmonics(q2, sph);
        std::cout << "n,l,m,center,Yp=" <<n<<" " << l<< " " <<m<< " (" << cart[0]<< " " << cart[1]<< " " << cart[2] <<") " << Yp << std::endl;
        std::cout << "n,l,-m,center,Yn=" <<n<<" " << l<< " " <<-1.0*m<< " (" << cart[0]<< " " << cart[1]<< " " << cart[2] <<") " << Yn << std::endl;
//        std::cout << "(Yp+(-1)^{m}Yn) /2=" <<n<<" " << l<< " " <<m<< " (" << cart[0]<< " " << cart[1]<< " " << cart[2] <<") " << sqrt(2)*(Yp+pow(-1,m)*Yn)/2.0 << std::endl;
//        std::cout << "n,l,m,center,Yp*=" <<n<<" " << l<< " " <<m<< " (" << cart[0]<< " " << cart[1]<< " " << cart[2] <<") " << conj(Yp) << std::endl;
//        std::cout << "n,l,m,center,-1^m (Yn)=" <<n<<" " << l<< " " <<m<< " (" << cart[0]<< " " << cart[1]<< " " << cart[2] <<") " << pow(-1,m)*Yn << std::endl;


    }
    return 0;
}