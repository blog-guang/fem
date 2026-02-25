#include "shape_tri3.h"

namespace fem {
namespace shape {

void Tri3ShapeFunction::evaluate(const Vec3& xi, Vector& N) const {
    Real xi_val = xi[0];
    Real eta = xi[1];
    
    N.resize(3);
    N[0] = 1.0 - xi_val - eta;  // N1
    N[1] = xi_val;               // N2
    N[2] = eta;                  // N3
}

void Tri3ShapeFunction::evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const {
    // 常数导数（线性单元）
    dN.resize(3, 2);
    
    // dN/dξ
    dN(0, 0) = -1.0;
    dN(1, 0) =  1.0;
    dN(2, 0) =  0.0;
    
    // dN/dη
    dN(0, 1) = -1.0;
    dN(1, 1) =  0.0;
    dN(2, 1) =  1.0;
}

void Tri3ShapeFunction::getGaussPoints(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    gaussTriangle(order, points, weights);
}

}  // namespace shape
}  // namespace fem
