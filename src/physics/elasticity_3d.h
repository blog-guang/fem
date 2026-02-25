/**
 * elasticity_3d.h - 3D 弹性力学物理模块
 * 
 * 使用形函数系统进行单元装配
 * 
 * 方程: σ = D ε, ε = ∇_s u
 * 
 * 使用方法:
 * ```cpp
 * Elasticity3D elast(E, nu);
 * Assembler assembler(model, 3);  // 矢量场 (u_x, u_y, u_z)
 * 
 * auto elem_func = [&](Index elem_id, const Mesh& mesh,
 *                      DenseMatrix& Ke, Vector& Fe) {
 *     elast.compute_element(elem_id, mesh, Ke, Fe);
 * };
 * 
 * assembler.assemble(elem_func);
 * ```
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "math/dense_matrix.h"
#include "math/vector.h"

namespace fem {
namespace physics {

/**
 * 3D 弹性力学物理模块
 * 
 * 支持单元类型:
 * - Tet4: 四面体4节点线性单元
 * - Brick8: 六面体8节点三线性单元
 */
class Elasticity3D {
public:
    /**
     * 构造函数
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     */
    Elasticity3D(Real youngs_modulus, Real poissons_ratio)
        : E_(youngs_modulus), nu_(poissons_ratio) {
        compute_D_matrix();
    }

    /**
     * 计算单元刚度矩阵和载荷向量
     * 
     * Ke_ij = ∫_Ω B^T D B dΩ
     * Fe = 0 (体力暂不考虑)
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param Ke 输出：单元刚度矩阵 (12x12 for Tet4, 24x24 for Brick8)
     * @param Fe 输出：单元载荷向量
     */
    void compute_element(Index elem_id, const Mesh& mesh,
                        DenseMatrix& Ke, Vector& Fe) const;

    /**
     * 仅计算刚度矩阵
     */
    void compute_stiffness(Index elem_id, const Mesh& mesh,
                          DenseMatrix& Ke) const;

    // ── 参数设置 ──
    void set_youngs_modulus(Real E) {
        E_ = E;
        compute_D_matrix();
    }

    void set_poissons_ratio(Real nu) {
        nu_ = nu;
        compute_D_matrix();
    }

    Real youngs_modulus() const { return E_; }
    Real poissons_ratio() const { return nu_; }

private:
    Real E_;              ///< 杨氏模量
    Real nu_;             ///< 泊松比

    DenseMatrix D_;       ///< 本构矩阵 (6x6)

    /**
     * 计算 3D 本构矩阵 D
     * 
     * D = E/((1+ν)(1-2ν)) * [[1-ν,  ν,    ν,    0,         0,         0        ],
     *                         [ν,    1-ν,  ν,    0,         0,         0        ],
     *                         [ν,    ν,    1-ν,  0,         0,         0        ],
     *                         [0,    0,    0,    (1-2ν)/2,  0,         0        ],
     *                         [0,    0,    0,    0,         (1-2ν)/2,  0        ],
     *                         [0,    0,    0,    0,         0,         (1-2ν)/2 ]]
     * 
     * 应力-应变关系:
     * {σ_xx, σ_yy, σ_zz, τ_yz, τ_xz, τ_xy}^T = D * {ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy}^T
     */
    void compute_D_matrix();
};

} // namespace physics
} // namespace fem
