#pragma once
#include "libslater.h"
#include <boost/math/special_functions.hpp>

namespace bm = boost::math;

namespace slater {


struct Spherical_Coordinates {
    double theta = 0;
    double phi = 0;
    double radius = 0;
    Spherical_Coordinates(const center_t &cartesian){
        {
            auto x = cartesian[0];
            auto y = cartesian[1];
            auto z = cartesian[2];
            auto pi = bm::constants::pi<double>();

            radius = sqrt(x*x + y*y + z*z);
            if (radius == 0){
                theta=0;
                phi =0;
            }
            else {
                theta = acos(z / radius);
                if (x > 0) {
                    phi = atan(y / x);
                } else if (x < 0) {
                    if (y >= 0) {
                        phi = atan(y / x) + pi;
                    } else {
                        phi = atan(y / x) - pi;
                    }
                }
                    // x==0
                else {
                    if (y > 0) {
                        phi = pi / 2;
                    } else if (y < 0) {
                        phi = -pi / 2;
                    }
                        //x==0, y==0, phi = undefined
                    else {
                        phi = 0;
                    }
                }
            }
        }
    };
    center_t  Cartesian_Coordinates() const
    {
        spatial_coordinate_t z = radius * cos(theta);
        spatial_coordinate_t y = radius * sin(theta) * sin(phi);
        spatial_coordinate_t x = radius * sin(theta) * cos(phi);
        center_t cartesian = {x,y,z};
        return cartesian;
    };
    Spherical_Coordinates() = default;
};

}