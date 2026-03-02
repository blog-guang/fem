/**
 * elasticity_nonlinear.cpp - 几何非线性弹性力学实现
 */

#include "physics/elasticity_nonlinear.h"
#include "shape/shape_function_factory.h"
#include "core/logger.h"
#include <stdexcept>

namespace fem {
namespace physics {

ElasticityNonlinear::ElasticityNonlinear(
    constitutive::Material* material,
    int dimension,
    bool enable_geometric_stiffness)
    : material_(material),
      dimension_(dimension),
      enable_geometric_stiffness_(enable_geometric_stiffness),
      stress_update_method_(StressUpdate::GREEN_LAGRANGE)
{
    if (!material_) {
        throw std::invalid_argument("ElasticityNonlinear: material cannot be null");
    }
    
    if (dimension_ != 2 && dimension_ != 3) {
        throw std::invalid_argument("ElasticityNonlinear: dimension must be 2 or 3");
    }
}

// ═══════════════════════════════════════════════════════════
// 几何非线性接口
// ═══════════════════════════════════════════════════════════

void ElasticityNonlinear::compute_element_nonlinear(
    Index elem_id,
    const Mesh& mesh,
    const Vector& u_current,
    DenseMatrix& Ke,
    Vector& Fe) const
{
    // 提取单元位移
    Vector u_elem = extract_element_displacement(elem_id, mesh, u_current);
    
    // 计算材料刚度矩阵 K_mat
    DenseMatrix K_mat;
    compute_material_stiffness(elem_id, mesh, u_elem, K_mat);
    
    // 计算几何刚度矩阵 K_geo（如果启用）
    if (enable_geometric_stiffness_) {
        DenseMatrix K_geo;
        compute_geometric_stiffness(elem_id, mesh, u_elem, K_geo);
        
        // K_t = K_mat + K_geo
        Ke = K_mat + K_geo;
    } else {
        Ke = K_mat;
    }
    
    // 计算内力向量
    Fe = compute_internal_force(elem_id, mesh, u_current);
}

Vector ElasticityNonlinear::compute_internal_force(
    Index elem_id,
    const Mesh& mesh,
    const Vector& u_current) const
{
    // TODO: 完整实现
    // 当前返回零向量（待实现）
    
    const auto& elem = mesh.element(elem_id);
    int num_nodes = elem.num_nodes();
    int elem_dofs = num_nodes * dimension_;
    
    Vector F_int(elem_dofs, 0.0);
    
    return F_int;
}

// ═══════════════════════════════════════════════════════════
// 线性接口（向后兼容）
// ═══════════════════════════════════════════════════════════

void ElasticityNonlinear::compute_element(
    Index elem_id,
    const Mesh& mesh,
    DenseMatrix& Ke,
    Vector& Fe) const
{
    // 小变形近似：使用零位移
    Vector u_zero(mesh.num_nodes() * dimension_, 0.0);
    
    compute_element_nonlinear(elem_id, mesh, u_zero, Ke, Fe);
}

// ═══════════════════════════════════════════════════════════
// 材料刚度矩阵
// ═══════════════════════════════════════════════════════════

void ElasticityNonlinear::compute_material_stiffness(
    Index elem_id,
    const Mesh& mesh,
    const Vector& u_elem,
    DenseMatrix& K_mat) const
{
    const auto& elem = mesh.element(elem_id);
    int num_nodes = elem.num_nodes();
    int elem_dofs = num_nodes * dimension_;
    
    // 获取节点坐标
    std::vector<Vec3> coords = get_element_coords(elem_id, mesh);
    
    // 获取形函数
    auto shape_func = shape::ShapeFunctionFactory::create(elem.type());
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = get_gauss_order(elem.type());
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    // 初始化刚度矩阵
    K_mat.resize(elem_dofs, elem_dofs);
    K_mat.fill(0.0);
    
    // 高斯积分
    for (std::size_t gp = 0; gp < gauss_points.size(); gp++) {
        const Vec3& xi = gauss_points[gp];
        Real w = weights[gp];
        
        // 计算形函数物理坐标导数
        DenseMatrix dN_dx;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dx);
        
        // 计算雅可比行列式
        DenseMatrix J_mat = shape_func->computeJacobian(xi, coords);
        Real det_J = compute_jacobian_determinant(J_mat);
        
        if (det_J <= 0.0) {
            throw std::runtime_error("Negative Jacobian detected in element " + 
                                   std::to_string(elem_id));
        }
        
        // 构造 B 矩阵
        DenseMatrix B;
        compute_B_matrix(dN_dx, num_nodes, B);
        
        // 获取材料切线刚度矩阵 D
        // TODO: 需要根据当前应变状态计算
        Vector strain(dimension_ == 3 ? 6 : 4, 0.0);  // 暂时使用零应变
        DenseMatrix D;
        
        constitutive::StateVariables state = material_->createState();
        material_->computeTangent(strain, D, state);
        
        // K_mat += w * det(J) * B^T * D * B
        DenseMatrix BT_D = B.transpose() * D;
        DenseMatrix BT_D_B = BT_D * B;
        
        K_mat = K_mat + BT_D_B * (w * det_J);
    }
}

// ═══════════════════════════════════════════════════════════
// 几何刚度矩阵
// ═══════════════════════════════════════════════════════════

void ElasticityNonlinear::compute_geometric_stiffness(
    Index elem_id,
    const Mesh& mesh,
    const Vector& u_elem,
    DenseMatrix& K_geo) const
{
    const auto& elem = mesh.element(elem_id);
    int num_nodes = elem.num_nodes();
    int elem_dofs = num_nodes * dimension_;
    
    // 获取节点坐标
    std::vector<Vec3> coords = get_element_coords(elem_id, mesh);
    
    // 获取形函数
    auto shape_func = shape::ShapeFunctionFactory::create(elem.type());
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = get_gauss_order(elem.type());
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    // 初始化几何刚度矩阵
    K_geo.resize(elem_dofs, elem_dofs);
    K_geo.fill(0.0);
    
    // 高斯积分
    for (std::size_t gp = 0; gp < gauss_points.size(); gp++) {
        const Vec3& xi = gauss_points[gp];
        Real w = weights[gp];
        
        // 计算形函数物理坐标导数
        DenseMatrix dN_dx;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dx);
        
        // 计算雅可比行列式
        DenseMatrix J_mat = shape_func->computeJacobian(xi, coords);
        Real det_J = compute_jacobian_determinant(J_mat);
        
        // 计算当前点的应力
        // TODO: 从应变计算应力
        DenseMatrix sigma(3, 3);
        sigma.fill(0.0);  // 暂时使用零应力
        
        // 构造 G 矩阵
        DenseMatrix G;
        compute_G_matrix(dN_dx, num_nodes, G);
        
        // 应力矩阵（Voigt 记号转换）
        // σ_matrix = [σ_xx  σ_xy  σ_xz]
        //            [σ_xy  σ_yy  σ_yz]
        //            [σ_xz  σ_yz  σ_zz]
        
        // K_geo += w * det(J) * G^T * σ * G
        // 简化实现：K_geo = ∫ Σ (σ_ij * ∂N_i/∂x_k * ∂N_j/∂x_k) dV
        
        for (int i = 0; i < num_nodes; i++) {
            for (int j = 0; j < num_nodes; j++) {
                for (int d = 0; d < dimension_; d++) {
                    Real k_geo_contrib = 0.0;
                    
                    // k_geo_ij += Σ_k (σ_dd * ∂N_i/∂x_k * ∂N_j/∂x_k)
                    for (int k = 0; k < dimension_; k++) {
                        k_geo_contrib += sigma(d, d) * dN_dx(i, k) * dN_dx(j, k);
                    }
                    
                    int I = i * dimension_ + d;
                    int J = j * dimension_ + d;
                    
                    K_geo(I, J) += w * det_J * k_geo_contrib;
                }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════
// B 矩阵和 G 矩阵
// ═══════════════════════════════════════════════════════════

void ElasticityNonlinear::compute_B_matrix(
    const DenseMatrix& dN_dx,
    int num_nodes,
    DenseMatrix& B) const
{
    int num_strain = (dimension_ == 3) ? 6 : 4;
    int elem_dofs = num_nodes * dimension_;
    
    B.resize(num_strain, elem_dofs);
    B.fill(0.0);
    
    if (dimension_ == 3) {
        // 3D B 矩阵
        for (int i = 0; i < num_nodes; i++) {
            int col = i * 3;
            
            // ε_xx = ∂u/∂x
            B(0, col + 0) = dN_dx(i, 0);
            
            // ε_yy = ∂v/∂y
            B(1, col + 1) = dN_dx(i, 1);
            
            // ε_zz = ∂w/∂z
            B(2, col + 2) = dN_dx(i, 2);
            
            // γ_xy = ∂u/∂y + ∂v/∂x
            B(3, col + 0) = dN_dx(i, 1);
            B(3, col + 1) = dN_dx(i, 0);
            
            // γ_yz = ∂v/∂z + ∂w/∂y
            B(4, col + 1) = dN_dx(i, 2);
            B(4, col + 2) = dN_dx(i, 1);
            
            // γ_xz = ∂u/∂z + ∂w/∂x
            B(5, col + 0) = dN_dx(i, 2);
            B(5, col + 2) = dN_dx(i, 0);
        }
    } else {
        // 2D B 矩阵（平面应变）
        for (int i = 0; i < num_nodes; i++) {
            int col = i * 2;
            
            // ε_xx = ∂u/∂x
            B(0, col + 0) = dN_dx(i, 0);
            
            // ε_yy = ∂v/∂y
            B(1, col + 1) = dN_dx(i, 1);
            
            // ε_zz = 0（平面应变）
            B(2, col + 0) = 0.0;
            B(2, col + 1) = 0.0;
            
            // γ_xy = ∂u/∂y + ∂v/∂x
            B(3, col + 0) = dN_dx(i, 1);
            B(3, col + 1) = dN_dx(i, 0);
        }
    }
}

void ElasticityNonlinear::compute_G_matrix(
    const DenseMatrix& dN_dx,
    int num_nodes,
    DenseMatrix& G) const
{
    int elem_dofs = num_nodes * dimension_;
    
    // G 矩阵大小：(dimension * dimension) x elem_dofs
    G.resize(dimension_ * dimension_, elem_dofs);
    G.fill(0.0);
    
    // G 矩阵定义：
    // ∇u = G * u_elem
    // 
    // ∂u/∂x  ∂u/∂y  ∂u/∂z
    // ∂v/∂x  ∂v/∂y  ∂v/∂z
    // ∂w/∂x  ∂w/∂y  ∂w/∂z
    
    for (int i = 0; i < num_nodes; i++) {
        for (int d = 0; d < dimension_; d++) {
            for (int k = 0; k < dimension_; k++) {
                int row = d * dimension_ + k;  // (d, k) 位置
                int col = i * dimension_ + d;  // 第 i 个节点的第 d 个 DOF
                
                G(row, col) = dN_dx(i, k);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

DenseMatrix ElasticityNonlinear::compute_element_deformation_gradient(
    const std::vector<Vec3>& coords,
    const Vector& u_elem,
    const DenseMatrix& dN_dxi,
    const std::vector<Real>& xi) const
{
    // TODO: 实现变形梯度计算
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    return F;
}

Vector ElasticityNonlinear::extract_element_displacement(
    Index elem_id,
    const Mesh& mesh,
    const Vector& u_global) const
{
    const auto& elem = mesh.element(elem_id);
    int num_nodes = elem.num_nodes();
    int elem_dofs = num_nodes * dimension_;
    
    Vector u_elem(elem_dofs);
    
    const auto& node_ids = elem.nodes();
    
    for (int i = 0; i < num_nodes; i++) {
        Index node_id = node_ids[i];
        
        for (int d = 0; d < dimension_; d++) {
            Index global_dof = node_id * dimension_ + d;
            
            if (global_dof < u_global.size()) {
                u_elem[i * dimension_ + d] = u_global[global_dof];
            } else {
                u_elem[i * dimension_ + d] = 0.0;
            }
        }
    }
    
    return u_elem;
}

}  // namespace physics
}  // namespace fem
