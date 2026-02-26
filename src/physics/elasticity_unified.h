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
 * // 方式 1: 简单接口（内部创建各向同性弹性材料）
 * ElasticityUnified elast_2d(E, nu, PlaneType::PlaneStress);
 * ElasticityUnified elast_3d(E, nu, true);  // 3D 模式
 * 
 * // 方式 2: 高级接口（使用自定义材料）
 * auto* material = new J2Plasticity(E, nu, sigma_y, H);  // 塑性材料
 * ElasticityUnified elast_custom(material, 3);  // 3D
 * 
 * // 装配
 * Assembler assembler(model, dim);
 * assembler.assemble([&](Index id, const Mesh& mesh, DenseMatrix& Ke, Vector& Fe) {
 *     elast.compute_element(id, mesh, Ke, Fe);
 * });
 * ```
 */

#pragma once

#include "physics/physics_base.h"
#include "material/material.h"
#include <memory>

namespace fem {
namespace physics {

/**
 * 平面问题类型 (2D)
 */
enum class PlaneType {
    PlaneStress,   ///< 平面应力 (薄板)
    PlaneStrain    ///< 平面应变 (厚板)
};

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
    // ═══════════════════════════════════════════════════════════
    // 构造函数 - 简单接口（内部创建各向同性弹性材料）
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 构造函数 (2D - 简单接口)
     * 内部创建 IsotropicElastic 材料
     * 
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     * @param plane_type 平面类型 (默认平面应力)
     */
    ElasticityUnified(Real youngs_modulus, Real poissons_ratio,
                     PlaneType plane_type = PlaneType::PlaneStress);
    
    /**
     * 构造函数 (3D - 简单接口)
     * 内部创建 IsotropicElastic 材料
     * 
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     * @param use_3d 设为 true 以启用 3D 模式
     */
    ElasticityUnified(Real youngs_modulus, Real poissons_ratio, bool use_3d);

    // ═══════════════════════════════════════════════════════════
    // 构造函数 - 高级接口（使用自定义材料）
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 构造函数 (高级接口 - 使用自定义材料)
     * 
     * @param material 材料本构模型指针（外部管理生命周期）
     * @param dimension 维度 (2 或 3)
     * 
     * 注意：
     * - material 指针由外部管理，ElasticityUnified 不拥有所有权
     * - material 必须在 ElasticityUnified 对象生命周期内有效
     * - 对于 2D 材料，确保 dimension 参数与材料构造时一致
     */
    ElasticityUnified(constitutive::Material* material, int dimension);
    
    /**
     * 析构函数
     * 如果内部创建了材料对象，会自动释放
     */
    ~ElasticityUnified();

    // ═══════════════════════════════════════════════════════════
    // 核心接口
    // ═══════════════════════════════════════════════════════════
    
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

    // ═══════════════════════════════════════════════════════════
    // 访问器
    // ═══════════════════════════════════════════════════════════
    
    int dimension() const { return dimension_; }
    bool is_2d() const { return dimension_ == 2; }
    const constitutive::Material* material() const { return material_; }
    
    /**
     * 获取材料参数（仅对简单构造函数有效）
     * 如果使用自定义材料，这些函数返回默认值或抛出警告
     */
    Real youngs_modulus() const;
    Real poissons_ratio() const;
    PlaneType plane_type() const;

private:
    constitutive::Material* material_;  ///< 材料本构模型指针
    bool own_material_;                 ///< 是否拥有材料对象（内部创建的需要释放）
    int dimension_;                     ///< 维度 (2 或 3)
    PlaneType plane_type_;              ///< 平面类型（仅 2D 简单构造函数有效）
};

} // namespace physics
} // namespace fem
