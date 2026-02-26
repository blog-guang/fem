/**
 * elasticity_unified.h - 统一的弹性力学物理模块
 * 
 * 方程: σ = D ε, ε = ∇_s u
 * 
 * 支持：
 * - 2D: 平面应力/平面应变 (Tri3, Quad4)
 * - 3D: 通用 3D 弹性 (Tet4, Brick8)
 * 
 * 特性：
 * - 集成材料本构模型系统 (Material)
 * - 无硬编码 D 矩阵，使用 Material::computeTangent()
 * - 支持弹性、塑性、损伤等任意材料模型
 * - 使用形函数系统进行通用装配
 * 
 * 使用方法:
 * ```cpp
 * // 创建材料
 * auto* material = new IsotropicElastic(E, nu, 2, true);  // 2D, plane_stress
 * // 或
 * auto* material = new J2Plasticity(E, nu, sigma_y, H, 3);  // 3D 塑性
 * 
 * // 创建 physics 对象
 * ElasticityUnified physics(material, 2);  // 2D
 * 
 * // 装配
 * Assembler assembler(model, dim);
 * assembler.assemble([&](Index id, const Mesh& mesh, DenseMatrix& Ke, Vector& Fe) {
 *     physics.compute_element(id, mesh, Ke, Fe);
 * });
 * ```
 */

#pragma once

#include "physics/physics_base.h"
#include "material/material.h"

namespace fem {
namespace physics {

/**
 * 统一的弹性力学物理模块
 * 
 * 特性：
 * - 自动识别 2D/3D
 * - 支持所有单元类型
 * - 集成材料本构模型系统
 * - 使用形函数系统
 * - 高斯积分精确计算
 * - 构造 B 矩阵，通过 Material 获取 D 矩阵
 */
class ElasticityUnified : public PhysicsBase {
public:
    /**
     * 构造函数（使用自定义材料）
     * 
     * @param material 材料本构模型指针（外部管理生命周期）
     * @param dimension 维度 (2 或 3)
     * 
     * 注意：
     * - material 指针由外部管理，ElasticityUnified 不拥有所有权
     * - material 必须在 ElasticityUnified 对象生命周期内有效
     * - 对于 2D 材料，确保 dimension 参数与材料构造时一致
     * 
     * 示例：
     * ```cpp
     * // 2D 平面应力
     * auto* mat = new IsotropicElastic(E, nu, 2, true);
     * ElasticityUnified physics(mat, 2);
     * 
     * // 3D 塑性
     * auto* mat = new J2Plasticity(E, nu, sigma_y, H, 3);
     * ElasticityUnified physics(mat, 3);
     * ```
     */
    ElasticityUnified(constitutive::Material* material, int dimension);
    
    /**
     * 析构函数
     */
    ~ElasticityUnified() = default;

    /**
     * 计算单元刚度矩阵和载荷向量
     * 
     * Ke = ∫_Ω B^T D B dΩ
     * Fe = 0 (体力暂不考虑)
     * 
     * D 矩阵通过 material_->computeTangent() 获取
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param Ke 输出：单元刚度矩阵
     * @param Fe 输出：单元载荷向量
     */
    void compute_element(Index elem_id, const Mesh& mesh,
                        DenseMatrix& Ke, Vector& Fe) const;

    /**
     * 仅计算刚度矩阵
     */
    void compute_stiffness(Index elem_id, const Mesh& mesh,
                          DenseMatrix& Ke) const;

    /**
     * 获取维度
     */
    int dimension() const { return dimension_; }
    
    /**
     * 是否为 2D
     */
    bool is_2d() const { return dimension_ == 2; }
    
    /**
     * 获取材料指针
     */
    const constitutive::Material* material() const { return material_; }

private:
    constitutive::Material* material_;  ///< 材料本构模型指针
    int dimension_;                     ///< 维度 (2 或 3)
};

} // namespace physics
} // namespace fem
