#include "physics/elasticity_3d.h"
#include "shape/shape_function_factory.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void Elasticity3D::compute_D_matrix() {
    D_ = DenseMatrix(6, 6, 0.0);
    
    Real factor = E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
    Real lambda = factor * nu_;         // 拉梅第一参数
    Real mu = factor * (1.0 - 2.0 * nu_) / 2.0;  // 剪切模量 (拉梅第二参数)
    
    // 对角块 (正应力)
    D_(0, 0) = factor * (1.0 - nu_);
    D_(1, 1) = factor * (1.0 - nu_);
    D_(2, 2) = factor * (1.0 - nu_);
    
    // 泊松耦合
    D_(0, 1) = lambda;
    D_(0, 2) = lambda;
    D_(1, 0) = lambda;
    D_(1, 2) = lambda;
    D_(2, 0) = lambda;
    D_(2, 1) = lambda;
    
    // 剪切块
    D_(3, 3) = mu;  // τ_yz
    D_(4, 4) = mu;  // τ_xz
    D_(5, 5) = mu;  // τ_xy
}

void Elasticity3D::compute_element(Index elem_id, const Mesh& mesh,
                                  DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    
    // 检查单元类型（支持 3D 单元）
    ElementType elem_type = elem.type();
    if (elem_type != ElementType::Tet4 && elem_type != ElementType::Brick8) {
        FEM_WARN("Elasticity3D: unsupported element type, skipping");
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
    Index n_dofs = n_nodes * 3;  // 3D: u_x, u_y, u_z per node
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
        DenseMatrix dN_dxyz;  // (n_nodes x 3)
        shape_func->computePhysicalDerivatives(xi, coords, dN_dxyz);
        
        // 计算雅可比行列式
        DenseMatrix J = shape_func->computeJacobian(xi, coords);
        Real det_J = J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1))
                   - J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0))
                   + J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));
        
        if (std::abs(det_J) < 1e-15) {
            FEM_WARN("Elasticity3D: degenerate element (det(J) near zero)");
            continue;
        }
        
        // 构造 B 矩阵 (6 x n_dofs)
        // ε = {ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy}^T
        // B = [[dN1/dx,  0,       0,       dN2/dx, ... ],
        //      [0,       dN1/dy,  0,       0,      ... ],
        //      [0,       0,       dN1/dz,  0,      ... ],
        //      [0,       dN1/dz,  dN1/dy,  0,      ... ],  // γ_yz = ∂u_y/∂z + ∂u_z/∂y
        //      [dN1/dz,  0,       dN1/dx,  ...     ... ],  // γ_xz = ∂u_x/∂z + ∂u_z/∂x
        //      [dN1/dy,  dN1/dx,  0,       ...     ... ]]  // γ_xy = ∂u_x/∂y + ∂u_y/∂x
        DenseMatrix B(6, n_dofs, 0.0);
        for (Index i = 0; i < n_nodes; ++i) {
            Index col_u = i * 3;      // u_x
            Index col_v = i * 3 + 1;  // u_y
            Index col_w = i * 3 + 2;  // u_z
            
            Real dN_dx = dN_dxyz(i, 0);
            Real dN_dy = dN_dxyz(i, 1);
            Real dN_dz = dN_dxyz(i, 2);
            
            // 正应变
            B(0, col_u) = dN_dx;  // ε_xx = ∂u_x/∂x
            B(1, col_v) = dN_dy;  // ε_yy = ∂u_y/∂y
            B(2, col_w) = dN_dz;  // ε_zz = ∂u_z/∂z
            
            // 剪应变
            B(3, col_v) = dN_dz;  // γ_yz = ∂u_y/∂z + ∂u_z/∂y
            B(3, col_w) = dN_dy;
            
            B(4, col_u) = dN_dz;  // γ_xz = ∂u_x/∂z + ∂u_z/∂x
            B(4, col_w) = dN_dx;
            
            B(5, col_u) = dN_dy;  // γ_xy = ∂u_x/∂y + ∂u_y/∂x
            B(5, col_v) = dN_dx;
        }
        
        // 计算刚度矩阵贡献: Ke += B^T * D * B * det(J) * w
        DenseMatrix DB = D_ * B;       // (6 x n_dofs)
        DenseMatrix Bt = B.transpose(); // (n_dofs x 6)
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

void Elasticity3D::compute_stiffness(Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke) const {
    Vector Fe;  // 占位
    compute_element(elem_id, mesh, Ke, Fe);
}

} // namespace physics
} // namespace fem
