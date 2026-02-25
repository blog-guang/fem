#include "shape_function.h"
#include <cmath>
#include <stdexcept>

namespace fem {
namespace shape {

DenseMatrix ShapeFunction::computeJacobian(
    const Vec3& xi,
    const std::vector<Vec3>& node_coords) const
{
    int dim = dimension();
    int nnodes = numNodes();
    
    if (node_coords.size() != static_cast<size_t>(nnodes)) {
        throw std::invalid_argument("node_coords size mismatch");
    }
    
    // 计算形函数导数
    DenseMatrix dN;
    evaluateDerivatives(xi, dN);
    
    // J = Σ (∂N_i/∂ξ_j) * x_i
    DenseMatrix J(dim, dim);
    J.zero();
    
    for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                J(j, k) += dN(i, k) * node_coords[i][j];
            }
        }
    }
    
    return J;
}

void ShapeFunction::computePhysicalDerivatives(
    const Vec3& xi,
    const std::vector<Vec3>& node_coords,
    DenseMatrix& dN_dx) const
{
    // 计算雅可比矩阵
    DenseMatrix J = computeJacobian(xi, node_coords);
    
    // 计算自然坐标系导数
    DenseMatrix dN_dxi;
    evaluateDerivatives(xi, dN_dxi);
    
    int dim = dimension();
    int nnodes = numNodes();
    
    // 计算 J^{-1}
    DenseMatrix J_inv(dim, dim);
    
    if (dim == 2) {
        // 2x2矩阵求逆
        Real det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
        if (std::abs(det) < 1e-15) {
            throw std::runtime_error("Singular Jacobian matrix");
        }
        J_inv(0, 0) =  J(1, 1) / det;
        J_inv(0, 1) = -J(0, 1) / det;
        J_inv(1, 0) = -J(1, 0) / det;
        J_inv(1, 1) =  J(0, 0) / det;
    } else if (dim == 3) {
        // 3x3矩阵求逆（伴随矩阵法）
        Real det = J(0, 0) * (J(1, 1)*J(2, 2) - J(1, 2)*J(2, 1))
                 - J(0, 1) * (J(1, 0)*J(2, 2) - J(1, 2)*J(2, 0))
                 + J(0, 2) * (J(1, 0)*J(2, 1) - J(1, 1)*J(2, 0));
        
        if (std::abs(det) < 1e-15) {
            throw std::runtime_error("Singular Jacobian matrix");
        }
        
        J_inv(0, 0) = (J(1, 1)*J(2, 2) - J(1, 2)*J(2, 1)) / det;
        J_inv(0, 1) = (J(0, 2)*J(2, 1) - J(0, 1)*J(2, 2)) / det;
        J_inv(0, 2) = (J(0, 1)*J(1, 2) - J(0, 2)*J(1, 1)) / det;
        
        J_inv(1, 0) = (J(1, 2)*J(2, 0) - J(1, 0)*J(2, 2)) / det;
        J_inv(1, 1) = (J(0, 0)*J(2, 2) - J(0, 2)*J(2, 0)) / det;
        J_inv(1, 2) = (J(0, 2)*J(1, 0) - J(0, 0)*J(1, 2)) / det;
        
        J_inv(2, 0) = (J(1, 0)*J(2, 1) - J(1, 1)*J(2, 0)) / det;
        J_inv(2, 1) = (J(0, 1)*J(2, 0) - J(0, 0)*J(2, 1)) / det;
        J_inv(2, 2) = (J(0, 0)*J(1, 1) - J(0, 1)*J(1, 0)) / det;
    }
    
    // dN/dx = dN/dξ * J^{-1}
    dN_dx.resize(nnodes, dim);
    dN_dx.zero();
    
    for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                dN_dx(i, j) += dN_dxi(i, k) * J_inv(k, j);
            }
        }
    }
}

Vec3 ShapeFunction::interpolate(
    const Vec3& xi,
    const std::vector<Vec3>& node_coords) const
{
    Vector N;
    evaluate(xi, N);
    
    int nnodes = numNodes();
    Vec3 x{0.0, 0.0, 0.0};
    
    for (int i = 0; i < nnodes; ++i) {
        x[0] += N[i] * node_coords[i][0];
        x[1] += N[i] * node_coords[i][1];
        x[2] += N[i] * node_coords[i][2];
    }
    
    return x;
}

// ═══ 2D Gauss积分点 ═══

void ShapeFunction2D::gaussLegendre2D(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    // 1D Gauss-Legendre点
    std::vector<Real> xi_1d, w_1d;
    
    if (order == 1) {
        xi_1d = {0.0};
        w_1d = {2.0};
    } else if (order == 2) {
        xi_1d = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
        w_1d = {1.0, 1.0};
    } else if (order == 3) {
        xi_1d = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
        w_1d = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else {
        throw std::invalid_argument("Unsupported Gauss order: " + std::to_string(order));
    }
    
    // 张量积生成2D点
    points.clear();
    weights.clear();
    
    for (size_t i = 0; i < xi_1d.size(); ++i) {
        for (size_t j = 0; j < xi_1d.size(); ++j) {
            points.push_back(Vec3{xi_1d[i], xi_1d[j], 0.0});
            weights.push_back(w_1d[i] * w_1d[j]);
        }
    }
}

void ShapeFunction2D::gaussTriangle(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    // 三角形高斯积分点（面积坐标）
    points.clear();
    weights.clear();
    
    if (order == 1) {
        // 1点积分：中心点
        points = {Vec3{1.0/3.0, 1.0/3.0, 0.0}};
        weights = {0.5};
    } else if (order == 2) {
        // 3点积分
        points = {
            Vec3{1.0/6.0, 1.0/6.0, 0.0},
            Vec3{2.0/3.0, 1.0/6.0, 0.0},
            Vec3{1.0/6.0, 2.0/3.0, 0.0}
        };
        weights = {1.0/6.0, 1.0/6.0, 1.0/6.0};
    } else if (order == 3) {
        // 6点积分
        Real a1 = 0.091576213509771;
        Real a2 = 0.445948490915965;
        Real w1 = 0.109951743655322;
        Real w2 = 0.223381589678011;
        
        points = {
            Vec3{a1, a1, 0.0},
            Vec3{1.0-2.0*a1, a1, 0.0},
            Vec3{a1, 1.0-2.0*a1, 0.0},
            Vec3{a2, a2, 0.0},
            Vec3{1.0-2.0*a2, a2, 0.0},
            Vec3{a2, 1.0-2.0*a2, 0.0}
        };
        weights = {w1, w1, w1, w2, w2, w2};
        
        // 归一化到参考三角形面积
        for (auto& w : weights) {
            w *= 0.5;
        }
    } else {
        throw std::invalid_argument("Unsupported triangle Gauss order: " + std::to_string(order));
    }
}

// ═══ 3D Gauss积分点 ═══

void ShapeFunction3D::gaussLegendre3D(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    // 1D Gauss-Legendre点
    std::vector<Real> xi_1d, w_1d;
    
    if (order == 1) {
        xi_1d = {0.0};
        w_1d = {2.0};
    } else if (order == 2) {
        xi_1d = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
        w_1d = {1.0, 1.0};
    } else if (order == 3) {
        xi_1d = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
        w_1d = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else {
        throw std::invalid_argument("Unsupported Gauss order: " + std::to_string(order));
    }
    
    // 张量积生成3D点
    points.clear();
    weights.clear();
    
    for (size_t i = 0; i < xi_1d.size(); ++i) {
        for (size_t j = 0; j < xi_1d.size(); ++j) {
            for (size_t k = 0; k < xi_1d.size(); ++k) {
                points.push_back(Vec3{xi_1d[i], xi_1d[j], xi_1d[k]});
                weights.push_back(w_1d[i] * w_1d[j] * w_1d[k]);
            }
        }
    }
}

void ShapeFunction3D::gaussTetrahedron(
    int order,
    std::vector<Vec3>& points,
    std::vector<Real>& weights) const
{
    // 四面体高斯积分点（体积坐标）
    points.clear();
    weights.clear();
    
    if (order == 1) {
        // 1点积分：中心点
        points = {Vec3{0.25, 0.25, 0.25}};
        weights = {1.0/6.0};
    } else if (order == 2) {
        // 4点积分
        Real a = 0.585410196624969;
        Real b = 0.138196601125011;
        
        points = {
            Vec3{a, b, b},
            Vec3{b, a, b},
            Vec3{b, b, a},
            Vec3{b, b, b}
        };
        weights = {0.25/6.0, 0.25/6.0, 0.25/6.0, 0.25/6.0};
    } else {
        throw std::invalid_argument("Unsupported tetrahedron Gauss order: " + std::to_string(order));
    }
}

}  // namespace shape
}  // namespace fem
