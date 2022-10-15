#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "slater-utils.h"
#include "libslater.h"

namespace bg = boost :: geometry ;

namespace slater
{

double distance(const center_t &A, const center_t &B)
{
    bg::model::point<double, 3, bg::cs::cartesian> A_point(A[0], A[1], A[2]);
    bg::model::point<double, 3, bg::cs::cartesian> B_point(B[0], B[1], B[2]);
    return bg::distance(A_point,B_point);
}


}
