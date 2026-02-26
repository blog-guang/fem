/**
 * heat_unified.h - 统一的热传导物理模块
 * 
 * 方程: -∇·(k∇T) = Q
 * 
 * 支持单元类型:
 * - 2D: Tri3, Quad4
 * - 3D: Tet4, Brick8
 * 
 * 使用形函数系统进行通用装配。
 * 
 * 使用方法:
 * ```cpp
 * HeatConductionUnified heat(k, Q);
 * Assembler assembler(model, 1);  // 标量场
 * 
 * auto elem_func = [&](Index elem_id, const Mesh& mesh,
 *                      DenseMatrix& Ke, Vector& Fe) {
 *     heat.compute_element(elem_id, mesh, Ke, Fe);
 * };
 * 
 * assembler.assemble(elem_func);
 * ```
 */

#pragma once

#include "physics/physics_base.h"

namespace fem {
namespace physics {

/**
 * 统一的热传导物理模块
 * 
 * 特性：
 * - 支持所有 2D/3D 单元
 * - 使用形函数系统（ShapeFunction）
 * - 高斯积分精确计算
 * - 无硬编码单元类型
 */
class HeatConductionUnified : public PhysicsBase {
public:
    /**
     * 构造函数
     * @param conductivity 导热系数 k (默认 1.0)
     * @param source 体热源 Q (默认 0.0)
     */
    HeatConductionUnified(Real conductivity = 1.0, Real source = 0.0)
        : k_(conductivity), Q_(source) {}

    /**
     * 计算单元刚度矩阵和载荷向量
     * 
     * Ke_ij = ∫_Ω k ∇Ni · ∇Nj dΩ
     * Fe_i  = ∫_Ω Q Ni dΩ
     * 
     * 使用高斯积分计算，支持所有单元类型。
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param Ke 输出：单元刚度矩阵 (n_nodes × n_nodes)
     * @param Fe 输出：单元载荷向量 (n_nodes)
     */
    void compute_element(Index elem_id, const Mesh& mesh,
                        DenseMatrix& Ke, Vector& Fe) const;

    /**
     * 仅计算刚度矩阵
     */
    void compute_stiffness(Index elem_id, const Mesh& mesh,
                          DenseMatrix& Ke) const;

    /**
     * 仅计算载荷向量
     */
    void compute_load(Index elem_id, const Mesh& mesh,
                     Vector& Fe) const;

    // ── 参数设置 ──
    void set_conductivity(Real k) { k_ = k; }
    void set_source(Real Q) { Q_ = Q; }

    Real conductivity() const { return k_; }
    Real source() const { return Q_; }

private:
    Real k_;  ///< 导热系数
    Real Q_;  ///< 体热源
};

} // namespace physics
} // namespace fem
