/**
 * elasticity_nonlinear.h - 几何非线性弹性力学模块
 * 
 * 理论：Updated Lagrangian 公式
 * 
 * 特性：
 * - 大变形分析
 * - 几何刚度矩阵 K_geo
 * - Green-Lagrange 应变
 * - 支持超弹性材料（Neo-Hookean, Mooney-Rivlin）
 * 
 * 切线刚度矩阵：
 *   K_t = K_mat + K_geo
 * 
 * 其中：
 *   K_mat = ∫ B_L^T D B_L dV  (材料刚度)
 *   K_geo = ∫ G^T σ G dV      (几何刚度)
 * 
 * 使用方法：
 * ```cpp
 * // 创建超弹性材料
 * auto mat = MaterialFactory::create("NeoHookean", {{"E", 10e6}, {"nu", 0.45}});
 * 
 * // 创建几何非线性物理模块
 * ElasticityNonlinear physics(mat.get(), 3);
 * 
 * // 在 Newton-Raphson 循环中
 * physics.compute_element(elem_id, mesh, u_current, Ke, Fe);
 * ```
 */

#pragma once

#include "physics/physics_base.h"
#include "material/material.h"
#include "core/kinematics.h"
#include "math/vector.h"
#include "math/dense_matrix.h"

namespace fem {
namespace physics {

/**
 * ElasticityNonlinear - 几何非线性弹性力学模块
 * 
 * 实现 Updated Lagrangian 公式：
 * 1. 从当前位移计算变形梯度 F
 * 2. 计算 Green-Lagrange 应变 E
 * 3. 通过材料本构获得 Cauchy 应力 σ
 * 4. 计算材料刚度矩阵 K_mat
 * 5. 计算几何刚度矩阵 K_geo
 * 6. 返回切线刚度 K_t = K_mat + K_geo
 * 
 * 注意：
 * - 需要当前位移向量 u_current
 * - 材料模型应使用大变形本构（如 Neo-Hookean）
 * - 小变形材料（如 IsotropicElastic）会退化为线性分析
 */
class ElasticityNonlinear : public PhysicsBase {
public:
    /**
     * 构造函数
     * 
     * @param material 材料本构模型（支持大变形）
     * @param dimension 维度（2 或 3）
     * @param enable_geometric_stiffness 是否启用几何刚度矩阵
     */
    ElasticityNonlinear(constitutive::Material* material, 
                       int dimension,
                       bool enable_geometric_stiffness = true);
    
    ~ElasticityNonlinear() = default;
    
    // ═══ 几何非线性接口 ═══
    
    /**
     * 计算单元刚度矩阵和内力向量（几何非线性版本）
     * 
     * 需要当前位移场：
     *   K_t = K_mat + K_geo
     *   F_int = ∫ B^T σ dV
     * 
     * @param elem_id 单元 ID
     * @param mesh 网格
     * @param u_current 当前位移向量（全局 DOFs）
     * @param Ke 输出：切线刚度矩阵 K_t
     * @param Fe 输出：内力向量 F_int
     */
    void compute_element_nonlinear(
        Index elem_id,
        const Mesh& mesh,
        const Vector& u_current,
        DenseMatrix& Ke,
        Vector& Fe
    ) const;
    
    /**
     * 计算单元内力向量（用于残差计算）
     * 
     * F_int = ∫ B^T σ dV
     * 
     * @param elem_id 单元 ID
     * @param mesh 网格
     * @param u_current 当前位移向量
     * @return F_int 内力向量
     */
    Vector compute_internal_force(
        Index elem_id,
        const Mesh& mesh,
        const Vector& u_current
    ) const;
    
    // ═══ 标准线性接口（向后兼容）═══
    
    /**
     * 线性版本（小变形近似）
     * 
     * K = ∫ B^T D B dV
     * Fe = 0
     * 
     * 注意：这是小变形近似，不使用当前位移
     */
    void compute_element(Index elem_id, const Mesh& mesh,
                        DenseMatrix& Ke, Vector& Fe) const override;
    
    // ═══ 设置 ═══
    
    /**
     * 启用/禁用几何刚度矩阵
     */
    void set_geometric_stiffness(bool enable) {
        enable_geometric_stiffness_ = enable;
    }
    
    bool is_geometric_stiffness_enabled() const {
        return enable_geometric_stiffness_;
    }
    
    /**
     * 设置应力更新方法
     */
    enum class StressUpdate {
        SMALL_STRAIN,      // 小应变：σ = D ε
        GREEN_LAGRANGE,    // Green-Lagrange 应变 + 超弹性
        TOTAL_LAGRANGIAN   // Total Lagrangian（未实现）
    };
    
    void set_stress_update_method(StressUpdate method) {
        stress_update_method_ = method;
    }
    
    // ═══ 访问器 ═══
    
    int dimension() const { return dimension_; }
    const constitutive::Material* material() const { return material_; }

private:
    // ═══ 数据成员 ═══
    
    constitutive::Material* material_;  ///< 材料本构模型
    int dimension_;                     ///< 维度（2/3）
    bool enable_geometric_stiffness_;   ///< 是否启用几何刚度
    StressUpdate stress_update_method_; ///< 应力更新方法
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 计算材料刚度矩阵
     * 
     * K_mat = ∫ B_L^T D B_L dV
     * 
     * B_L 是线性应变-位移矩阵
     */
    void compute_material_stiffness(
        Index elem_id,
        const Mesh& mesh,
        const Vector& u_elem,
        DenseMatrix& K_mat
    ) const;
    
    /**
     * 计算几何刚度矩阵
     * 
     * K_geo = ∫ G^T σ G dV
     * 
     * G 是位移梯度-位移矩阵
     * σ 是 Cauchy 应力
     */
    void compute_geometric_stiffness(
        Index elem_id,
        const Mesh& mesh,
        const Vector& u_elem,
        DenseMatrix& K_geo
    ) const;
    
    /**
     * 从单元位移计算单元变形梯度
     * 
     * @param coords 单元节点坐标
     * @param u_elem 单元位移
     * @param dN_dxi 形函数导数（参考坐标）
     * @param xi 高斯点位置
     * @return F 变形梯度（3x3）
     */
    DenseMatrix compute_element_deformation_gradient(
        const std::vector<Vec3>& coords,
        const Vector& u_elem,
        const DenseMatrix& dN_dxi,
        const std::vector<Real>& xi
    ) const;
    
    /**
     * 计算线性应变-位移矩阵 B_L
     * 
     * ε = B_L * u_elem
     */
    void compute_B_matrix(
        const DenseMatrix& dN_dx,
        int num_nodes,
        DenseMatrix& B
    ) const;
    
    /**
     * 计算位移梯度矩阵 G
     * 
     * ∇u = G * u_elem
     */
    void compute_G_matrix(
        const DenseMatrix& dN_dx,
        int num_nodes,
        DenseMatrix& G
    ) const;
    
    /**
     * 提取单元位移向量
     * 
     * @param elem_id 单元 ID
     * @param mesh 网格
     * @param u_global 全局位移向量
     * @return u_elem 单元位移向量
     */
    Vector extract_element_displacement(
        Index elem_id,
        const Mesh& mesh,
        const Vector& u_global
    ) const;
};

}  // namespace physics
}  // namespace fem
