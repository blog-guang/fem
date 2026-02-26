#include "physics/elasticity_unified.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

// ═══════════════════════════════════════════════════════════
// 构造函数实现
// ═══════════════════════════════════════════════════════════

ElasticityUnified::ElasticityUnified(constitutive::Material* material, int dimension)
    : material_(material), dimension_(dimension) {
    if (!material_) {
        FEM_ERROR("ElasticityUnified: material pointer is null");
    }
    if (dimension_ != 2 && dimension_ != 3) {
        FEM_ERROR("ElasticityUnified: dimension must be 2 or 3");
    }
}

// ═══════════════════════════════════════════════════════════
// 核心计算函数
// ═══════════════════════════════════════════════════════════

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
    int elem_dim = shape_func->dimension();
    
    if (dimension_ != elem_dim) {
        FEM_WARN("ElasticityUnified: dimension mismatch (physics: " 
                 + std::to_string(dimension_) + ", element: " 
                 + std::to_string(elem_dim) + ")");
        return;
    }
    
    int dofs_per_node = dimension_;
    
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
    
    // 创建状态变量（弹性材料为空，塑性材料需要）
    constitutive::StateVariables state = material_->createState();
    
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
        if (dimension_ == 2) {
            B = build_B_matrix_2D(dN_dxyz);
        } else {
            B = build_B_matrix_3D(dN_dxyz);
        }
        
        // 获取材料切线刚度矩阵 D
        // 注意：对于线性弹性，D 不依赖应变，但对于非线性材料可能需要当前应变状态
        Vector strain_dummy;  // 弹性材料不需要应变信息
        DenseMatrix D;
        material_->computeTangent(strain_dummy, D, state);
        
        // 装配刚度矩阵: Ke += B^T * D * B * dV
        DenseMatrix DB = D * B;  // D * B (strain_size × n_dofs)
        
        for (Index i = 0; i < n_dofs; ++i) {
            for (Index j = 0; j < n_dofs; ++j) {
                Real sum = 0.0;
                Index strain_size = D.rows();
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
