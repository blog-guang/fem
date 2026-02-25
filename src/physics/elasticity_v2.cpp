#include "physics/elasticity_v2.h"
#include "shape/shape_function_factory.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void Elasticity2D::compute_D_matrix() {
    D_ = DenseMatrix(3, 3, 0.0);
    
    if (plane_type_ == PlaneType::PlaneStress) {
        // 平面应力
        Real factor = E_ / (1.0 - nu_ * nu_);
        D_(0, 0) = factor;
        D_(0, 1) = factor * nu_;
        D_(1, 0) = factor * nu_;
        D_(1, 1) = factor;
        D_(2, 2) = factor * (1.0 - nu_) / 2.0;
    } else {
        // 平面应变
        Real factor = E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
        D_(0, 0) = factor * (1.0 - nu_);
        D_(0, 1) = factor * nu_;
        D_(1, 0) = factor * nu_;
        D_(1, 1) = factor * (1.0 - nu_);
        D_(2, 2) = factor * (1.0 - 2.0 * nu_) / 2.0;
    }
}

void Elasticity2D::compute_element(Index elem_id, const Mesh& mesh,
                                  DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    
    // 检查单元类型（支持 2D 单元）
    ElementType elem_type = elem.type();
    if (elem_type != ElementType::Tri3 && elem_type != ElementType::Quad4) {
        FEM_WARN("Elasticity2D: unsupported element type, skipping");
        return;
    }
    
    // 获取单元节点坐标
    const auto& nodes = elem.nodes();
    Index n_nodes = nodes.size();
    std::vector<Vec3> coords(n_nodes);
    for (Index i = 0; i < n_nodes; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    // 创建形函数对象
    auto shape_func = shape::ShapeFunctionFactory::create(elem_type);
    
    // 获取高斯积分点（阶数2：精确积分线性应变）
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    shape_func->getGaussPoints(2, gauss_points, weights);
    
    // 初始化单元矩阵
    Index n_dofs = n_nodes * 2;  // 2D: u_x, u_y per node
    Ke.resize(n_dofs, n_dofs);
    Ke.zero();
    Fe.resize(n_dofs);
    for (Index i = 0; i < n_dofs; ++i) {
        Fe[i] = 0.0;
    }
    
    // 高斯积分循环
    for (size_t gp = 0; gp < gauss_points.size(); ++gp) {
        const Vec3& xi = gauss_points[gp];
        Real w = weights[gp];
        
        // 计算形函数导数（物理坐标系）
        DenseMatrix dN_dx;  // (n_nodes x 2)
        shape_func->computePhysicalDerivatives(xi, coords, dN_dx);
        
        // 计算雅可比行列式
        DenseMatrix J = shape_func->computeJacobian(xi, coords);
        Real det_J = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
        
        if (std::abs(det_J) < 1e-15) {
            FEM_WARN("Elasticity2D: degenerate element (det(J) near zero)");
            continue;
        }
        
        // 构造 B 矩阵 (3 x n_dofs)
        // B = [[dN1/dx,  0,       dN2/dx, 0,       ... ],
        //      [0,       dN1/dy,  0,      dN2/dy,  ... ],
        //      [dN1/dy,  dN1/dx,  dN2/dy, dN2/dx,  ... ]]
        DenseMatrix B(3, n_dofs, 0.0);
        for (Index i = 0; i < n_nodes; ++i) {
            Index col_u = i * 2;      // u_x 列
            Index col_v = i * 2 + 1;  // u_y 列
            
            B(0, col_u) = dN_dx(i, 0);  // ε_xx = ∂u_x/∂x
            B(1, col_v) = dN_dx(i, 1);  // ε_yy = ∂u_y/∂y
            B(2, col_u) = dN_dx(i, 1);  // γ_xy = ∂u_x/∂y + ∂u_y/∂x
            B(2, col_v) = dN_dx(i, 0);
        }
        
        // 计算刚度矩阵贡献: Ke += B^T * D * B * det(J) * w
        DenseMatrix DB = D_ * B;       // (3 x n_dofs)
        DenseMatrix Bt = B.transpose(); // (n_dofs x 3)
        DenseMatrix BtDB = Bt * DB;     // (n_dofs x n_dofs)
        
        Real factor = det_J * w;
        for (Index i = 0; i < n_dofs; ++i) {
            for (Index j = 0; j < n_dofs; ++j) {
                Ke(i, j) += BtDB(i, j) * factor;
            }
        }
        
        // 体力载荷 (暂不实现)
        // Fe += N^T * f * det(J) * w
    }
}

void Elasticity2D::compute_stiffness(Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke) const {
    Vector Fe;  // 占位
    compute_element(elem_id, mesh, Ke, Fe);
}

Real Elasticity2D::compute_tri3_B_matrix(const Vec3* coords, DenseMatrix& B) const {
    // 保留向后兼容（已弃用）
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];
    
    // 计算面积
    Real detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    Real area = 0.5 * std::abs(detJ);
    
    if (area < 1e-15) {
        return area;
    }
    
    // 计算形函数梯度
    Real dN_dx[3], dN_dy[3];
    
    // ∇N0
    dN_dx[0] = (y1 - y2) / (2.0 * area);
    dN_dy[0] = (x2 - x1) / (2.0 * area);
    
    // ∇N1
    dN_dx[1] = (y2 - y0) / (2.0 * area);
    dN_dy[1] = (x0 - x2) / (2.0 * area);
    
    // ∇N2
    dN_dx[2] = (y0 - y1) / (2.0 * area);
    dN_dy[2] = (x1 - x0) / (2.0 * area);
    
    // 构造 B 矩阵 (3x6)
    B.zero();
    
    for (int i = 0; i < 3; ++i) {
        int col_u = i * 2;      // u_x 列
        int col_v = i * 2 + 1;  // u_y 列
        
        B(0, col_u) = dN_dx[i];  // ε_xx = ∂u_x/∂x
        B(1, col_v) = dN_dy[i];  // ε_yy = ∂u_y/∂y
        B(2, col_u) = dN_dy[i];  // γ_xy = ∂u_x/∂y + ∂u_y/∂x
        B(2, col_v) = dN_dx[i];
    }
    
    return area;
}

} // namespace physics
} // namespace fem
