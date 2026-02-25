#include "shape_brick8.h"

namespace fem {
namespace shape {

void Brick8ShapeFunction::evaluate(const Vec3& xi, Vector& N) const {
    Real xi_val = xi[0], eta = xi[1], zeta = xi[2];
    
    N.resize(8);
    
    // 节点自然坐标：(-1,-1,-1), (+1,-1,-1), ...
    Real xi_nodes[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    Real eta_nodes[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    Real zeta_nodes[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    
    for (int i = 0; i < 8; ++i) {
        N[i] = 0.125 * (1.0 + xi_val * xi_nodes[i]) *
                       (1.0 + eta * eta_nodes[i]) *
                       (1.0 + zeta * zeta_nodes[i]);
    }
}

void Brick8ShapeFunction::evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const {
    Real xi_val = xi[0], eta = xi[1], zeta = xi[2];
    
    dN.resize(8, 3);
    
    Real xi_nodes[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    Real eta_nodes[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    Real zeta_nodes[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    
    for (int i = 0; i < 8; ++i) {
        // dN/dξ
        dN(i, 0) = 0.125 * xi_nodes[i] *
                   (1.0 + eta * eta_nodes[i]) *
                   (1.0 + zeta * zeta_nodes[i]);
        
        // dN/dη
        dN(i, 1) = 0.125 * (1.0 + xi_val * xi_nodes[i]) *
                   eta_nodes[i] *
                   (1.0 + zeta * zeta_nodes[i]);
        
        // dN/dζ
        dN(i, 2) = 0.125 * (1.0 + xi_val * xi_nodes[i]) *
                   (1.0 + eta * eta_nodes[i]) *
                   zeta_nodes[i];
    }
}

void Brick8ShapeFunction::getGaussPoints(int order, std::vector<Vec3>& points,
                                        std::vector<Real>& weights) const {
    gaussLegendre3D(order, points, weights);
}

}  // namespace shape
}  // namespace fem
