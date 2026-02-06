/**
 * heat.h - 热传导物理模块
 * 
 * 方程: -∇·(k∇u) = Q
 * 
 * k: 导热系数
 * Q: 体热源
 * u: 温度场 (标量)
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "math/dense_matrix.h"
#include "math/vector.h"

namespace fem {
namespace physics {

/**
 * 热传导物理模块
 * 
 * 使用方法:
 * ```cpp
 * HeatConduction heat(k, Q);
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
class HeatConduction {
public:
    /**
     * 构造函数
     * @param conductivity 导热系数 k (默认 1.0)
     * @param source 体热源 Q (默认 0.0)
     */
    HeatConduction(Real conductivity = 1.0, Real source = 0.0)
        : k_(conductivity), Q_(source) {}

    /**
     * 计算单元刚度矩阵和载荷向量
     * 
     * Ke_ij = ∫_Ω k ∇Ni · ∇Nj dΩ
     * Fe_i = ∫_Ω Q Ni dΩ
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

    /**
     * 计算三角形单元的形函数梯度
     * @param coords 节点坐标 [3]
     * @param grad 输出：梯度 [3][2] (dN/dx, dN/dy)
     * @return 单元面积
     */
    Real compute_tri3_gradients(const Vec3* coords, Real grad[][2]) const;
};

} // namespace physics
} // namespace fem
