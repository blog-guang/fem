/**
 * elasticity_v2.h - 弹性力学物理模块 (新架构)
 * 
 * 2D 平面应力/平面应变问题
 * 
 * 方程: σ = D ε, ε = ∇_s u
 * 
 * 使用方法:
 * ```cpp
 * Elasticity2D elast(E, nu, PlaneStress);
 * Assembler assembler(model, 2);  // 矢量场 (u_x, u_y)
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
 * 平面问题类型
 */
enum class PlaneType {
    PlaneStress,   ///< 平面应力 (薄板)
    PlaneStrain    ///< 平面应变 (厚板)
};

/**
 * 2D 弹性力学物理模块
 */
class Elasticity2D {
public:
    /**
     * 构造函数
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     * @param plane_type 平面类型 (平面应力/平面应变)
     */
    Elasticity2D(Real youngs_modulus, Real poissons_ratio,
                PlaneType plane_type = PlaneType::PlaneStress)
        : E_(youngs_modulus), nu_(poissons_ratio), plane_type_(plane_type) {
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
     * @param Ke 输出：单元刚度矩阵 (6x6 for Tri3)
     * @param Fe 输出：单元载荷向量 (6 for Tri3)
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

    void set_plane_type(PlaneType type) {
        plane_type_ = type;
        compute_D_matrix();
    }

    Real youngs_modulus() const { return E_; }
    Real poissons_ratio() const { return nu_; }
    PlaneType plane_type() const { return plane_type_; }

private:
    Real E_;              ///< 杨氏模量
    Real nu_;             ///< 泊松比
    PlaneType plane_type_; ///< 平面类型

    DenseMatrix D_;       ///< 本构矩阵 (3x3)

    /**
     * 计算本构矩阵 D
     * 
     * 平面应力:
     * D = E/(1-ν²) * [[1,   ν,   0         ],
     *                  [ν,   1,   0         ],
     *                  [0,   0,   (1-ν)/2 ]]
     * 
     * 平面应变:
     * D = E/((1+ν)(1-2ν)) * [[1-ν,  ν,     0           ],
     *                         [ν,    1-ν,   0           ],
     *                         [0,    0,     (1-2ν)/2  ]]
     */
    void compute_D_matrix();

    /**
     * 计算三角形单元的 B 矩阵 (应变-位移矩阵)
     * 
     * B = [[dN1/dx,  0,       dN2/dx, 0,       dN3/dx, 0      ],
     *      [0,       dN1/dy,  0,      dN2/dy,  0,      dN3/dy ],
     *      [dN1/dy,  dN1/dx,  dN2/dy, dN2/dx,  dN3/dy, dN3/dx]]
     * 
     * @param coords 节点坐标 [3]
     * @param B 输出：B 矩阵 (3x6)
     * @return 单元面积
     */
    Real compute_tri3_B_matrix(const Vec3* coords, DenseMatrix& B) const;
};

} // namespace physics
} // namespace fem
