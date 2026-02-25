#include "shape_tet4.h"

namespace fem {
namespace shape {

void Tet4ShapeFunction::evaluate(const Vec3& xi, Vector& N) const {
    Real xi_val = xi[0], eta = xi[1], zeta = xi[2];
    
    N.resize(4);
    N[0] = 1.0 - xi_val - eta - zeta;
    N[1] = xi_val;
    N[2] = eta;
    N[3] = zeta;
}

void Tet4ShapeFunction::evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const {
    dN.resize(4, 3);
    
    // dN/dξ
    dN(0, 0) = -1.0;  dN(0, 1) = -1.0;  dN(0, 2) = -1.0;
    dN(1, 0) =  1.0;  dN(1, 1) =  0.0;  dN(1, 2) =  0.0;
    dN(2, 0) =  0.0;  dN(2, 1) =  1.0;  dN(2, 2) =  0.0;
    dN(3, 0) =  0.0;  dN(3, 1) =  0.0;  dN(3, 2) =  1.0;
}

void Tet4ShapeFunction::getGaussPoints(int order, std::vector<Vec3>& points,
                                       std::vector<Real>& weights) const {
    gaussTetrahedron(order, points, weights);
}

}  // namespace shape
}  // namespace fem
