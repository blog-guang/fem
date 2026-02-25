#include "shape_quad4.h"

namespace fem {
namespace shape {

void Quad4ShapeFunction::evaluate(const Vec3& xi, Vector& N) const {
    Real xi_val = xi[0];
    Real eta = xi[1];
    
    N.resize(4);
    
    // 双线性形函数
    N[0] = 0.25 * (1.0 - xi_val) * (1.0 - eta);  // 节点0: (-1,-1)
    N[1] = 0.25 * (1.0 + xi_val) * (1.0 - eta);  // 节点1: (+1,-1)
    N[2] = 0.25 * (1.0 + xi_val) * (1.0 + eta);  // 节点2: (+1,+1)
    N[3] = 0.25 * (1.0 - xi_val) * (1.0 + eta);  // 节点3: (-1,+1)
}

void Quad4ShapeFunction::evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const {
    Real xi_val = xi[0];
    Real eta = xi[1];
    
    dN.resize(4, 2);
    
    // dN/dξ
    dN(0, 0) = -0.25 * (1.0 - eta);
    dN(1, 0) =  0.25 * (1.0 - eta);
    dN(2, 0) =  0.25 * (1.0 + eta);
    dN(3, 0) = -0.25 * (1.0 + eta);
    
    // dN/dη
    dN(0, 1) = -0.25 * (1.0 - xi_val);
    dN(1, 1) = -0.25 * (1.0 + xi_val);
    dN(2, 1) =  0.25 * (1.0 + xi_val);
    dN(3, 1) =  0.25 * (1.0 - xi_val);
}

void Quad4ShapeFunction::getGaussPoints(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    gaussLegendre2D(order, points, weights);
}

}  // namespace shape
}  // namespace fem
