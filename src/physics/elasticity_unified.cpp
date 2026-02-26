#include "physics/elasticity_unified.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void ElasticityUnified::compute_D_matrix_2D() {
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

void ElasticityUnified::compute_D_matrix_3D() {
    D_ = DenseMatrix(6, 6, 0.0);
    
    Real factor = E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
    Real lambda = factor * nu_;
    Real mu = factor * (1.0 - 2.0 * nu_) / 2.0;
    
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

void ElasticityUnified::compute_element(Index elem_id, const Mesh& mesh,
                                       DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    ElementType elem_type = elem.type();
    
    // 获取节点坐标
    std::vector<Vec3> coords = get_element_coords(elem_id, mesh);
    Index n_nodes = coords.size();
    
    // 创建形函数对象
    auto shape_func = shape::ShapeFunctionFactory::create(elem_type);
    if (!shape_func) {
        FEM_WARN("ElasticityUnified: unsupported element type");
        return;
    }
    
    // 检查维度一致性
    int dim = shape_func->dimension();
    int dofs_per_node = dim;
    
    if (is_2d_ && dim != 2) {
        FEM_WARN("ElasticityUnified: 2D mode but element is not 2D");
        return;
    }
    if (!is_2d_ && dim != 3) {
        FEM_WARN("ElasticityUnified: 3D mode but element is not 3D");
        return;
    }
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = get_gauss_order(elem_type);
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    // 初始化单元矩阵
    Index n_dofs = n_nodes * dofs_per_node;
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
        
        // 计算形函数物理坐标导数
        DenseMatrix dN_dxyz;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dxyz);
        
        // 计算雅可比行列式
        DenseMatrix J = shape_func->computeJacobian(xi, coords);
        Real det_J = compute_jacobian_determinant(J);
        
        if (std::abs(det_J) < 1e-15) {
            FEM_WARN("ElasticityUnified: degenerate element (det(J) near zero)");
            continue;
        }
        
        Real dV = w * std::abs(det_J);
        
        // 构造 B 矩阵
        DenseMatrix B;
        if (dim == 2) {
            B = build_B_matrix_2D(dN_dxyz);
        } else {
            B = build_B_matrix_3D(dN_dxyz);
        }
        
        // 装配刚度矩阵: Ke += B^T * D * B * dV
        // K_ij = Σ_gp B^T[strain_comp, i] * D[strain_comp, strain_comp2] * B[strain_comp2, j] * dV
        
        DenseMatrix DB = D_ * B;  // D * B (strain_size × n_dofs)
        
        for (Index i = 0; i < n_dofs; ++i) {
            for (Index j = 0; j < n_dofs; ++j) {
                Real sum = 0.0;
                Index strain_size = D_.rows();
                for (Index k = 0; k < strain_size; ++k) {
                    sum += B(k, i) * DB(k, j);
                }
                Ke(i, j) += sum * dV;
            }
        }
        
        // Fe = 0 (体力暂不考虑)
    }
}

void ElasticityUnified::compute_stiffness(Index elem_id, const Mesh& mesh,
                                         DenseMatrix& Ke) const {
    Vector Fe_dummy;
    compute_element(elem_id, mesh, Ke, Fe_dummy);
}

} // namespace physics
} // namespace fem
