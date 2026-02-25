/**
 * elasticity_material.h - 基于材料本构的弹性力学模块
 * 
 * 集成 Material 类的弹性力学求解器
 * 支持弹性和塑性材料
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "math/dense_matrix.h"
#include "math/vector.h"
#include "material/material.h"
#include <memory>
#include <vector>

namespace fem {
namespace physics {

using constitutive::Material;
using constitutive::StateVariables;

/**
 * 基于材料本构的2D弹性力学模块
 * 
 * 特性：
 * - 使用 Material 抽象接口
 * - 支持弹性和塑性材料
 * - 管理每个积分点的状态变量
 * - 增量法求解（支持非线性）
 */
class ElasticityWithMaterial {
public:
    /**
     * 构造函数
     * @param material 材料本构模型（共享指针）
     * @param thickness 厚度（2D问题）
     */
    ElasticityWithMaterial(std::shared_ptr<Material> material, 
                          Real thickness = 1.0)
        : material_(material), thickness_(thickness) {}
    
    // ═══ 初始化 ═══
    
    /**
     * 初始化状态变量（在装配前调用）
     * @param num_elements 单元总数
     */
    void initialize(std::size_t num_elements);
    
    // ═══ 单元计算 ═══
    
    /**
     * 计算单元刚度矩阵（切线刚度）
     * 
     * Ke = ∫_Ω B^T D B dΩ
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param Ke 输出：单元刚度矩阵
     */
    void compute_stiffness(Index elem_id, const Mesh& mesh,
                          DenseMatrix& Ke) const;
    
    /**
     * 更新应力（给定位移增量）
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param u_inc 位移增量向量（全局）
     * @param Fe 输出：单元内力向量
     */
    void update_stress(Index elem_id, const Mesh& mesh,
                      const Vector& u_inc, Vector& Fe);
    
    /**
     * 计算单元刚度和内力（用于Newton-Raphson）
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param u_total 总位移向量（全局）
     * @param Ke 输出：切线刚度矩阵
     * @param Fe 输出：内力向量
     */
    void compute_element(Index elem_id, const Mesh& mesh,
                        const Vector& u_total,
                        DenseMatrix& Ke, Vector& Fe);
    
    // ═══ 后处理 ═══
    
    /**
     * 计算单元应力（在积分点）
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param u_total 总位移向量
     * @return 应力向量（Voigt记号）
     */
    Vector compute_stress(Index elem_id, const Mesh& mesh,
                         const Vector& u_total) const;
    
    /**
     * 获取单元状态变量（只读）
     */
    const StateVariables& get_state(Index elem_id) const;
    
    /**
     * 重置所有状态变量
     */
    void reset();
    
    // ═══ 参数访问 ═══
    
    std::shared_ptr<Material> material() const { return material_; }
    Real thickness() const { return thickness_; }
    
    void set_thickness(Real t) { thickness_ = t; }

private:
    std::shared_ptr<Material> material_;  ///< 材料本构
    Real thickness_;                      ///< 厚度
    
    // 每个单元的状态变量（一个积分点）
    mutable std::vector<StateVariables> element_states_;
    mutable std::vector<Vector> element_stresses_;  // 当前应力
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 计算 Tri3 单元的 B 矩阵
     * 
     * B 矩阵：应变-位移关系
     * ε = B * u_e
     * 
     * @param coords 三个节点坐标
     * @param B 输出：B矩阵 (3x6)
     * @return 单元面积
     */
    Real compute_B_matrix(const Vec3 coords[3], DenseMatrix& B) const;
};

}  // namespace physics
}  // namespace fem
