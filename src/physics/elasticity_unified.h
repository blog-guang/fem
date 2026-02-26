/**
 * elasticity_unified.h - 统一的弹性力学物理模块
 * 
 * 方程: σ = D ε, ε = ∇_s u
 * 
 * 支持：
 * - 2D: 平面应力/平面应变 (Tri3, Quad4)
 * - 3D: 通用 3D 弹性 (Tet4, Brick8)
 * 
 * 使用形函数系统进行通用装配。
 * 
 * 使用方法:
 * ```cpp
 * // 2D
 * ElasticityUnified elast_2d(E, nu, PlaneType::PlaneStress);
 * Assembler assembler(model, 2);  // 矢量场 (u_x, u_y)
 * 
 * // 3D
 * ElasticityUnified elast_3d(E, nu);
 * Assembler assembler(model, 3);  // 矢量场 (u_x, u_y, u_z)
 * 
 * assembler.assemble([&](Index id, const Mesh& mesh, DenseMatrix& Ke, Vector& Fe) {
 *     elast.compute_element(id, mesh, Ke, Fe);
 * });
 * ```
 */

#pragma once

#include "physics/physics_base.h"

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
 * - 使用形函数系统
 * - 高斯积分精确计算
 * - 构造 B 矩阵和本构矩阵 D
 */
class ElasticityUnified : public PhysicsBase {
public:
    /**
     * 构造函数 (2D)
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     * @param plane_type 平面类型 (默认平面应力)
     */
    ElasticityUnified(Real youngs_modulus, Real poissons_ratio,
                     PlaneType plane_type = PlaneType::PlaneStress)
        : E_(youngs_modulus), nu_(poissons_ratio), 
          is_2d_(true), plane_type_(plane_type) {
        compute_D_matrix_2D();
    }
    
    /**
     * 构造函数 (3D - 通过标签区分)
     * @param youngs_modulus 杨氏模量 E
     * @param poissons_ratio 泊松比 ν
     * @param use_3d 设为 true 以启用 3D 模式
     */
    ElasticityUnified(Real youngs_modulus, Real poissons_ratio, bool use_3d)
        : E_(youngs_modulus), nu_(poissons_ratio), is_2d_(!use_3d) {
        if (use_3d) {
            compute_D_matrix_3D();
        } else {
            plane_type_ = PlaneType::PlaneStress;
            compute_D_matrix_2D();
        }
    }

    /**
     * 计算单元刚度矩阵和载荷向量
     * 
     * Ke = ∫_Ω B^T D B dΩ
     * Fe = 0 (体力暂不考虑)
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

    // ── 参数设置 ──
    void set_youngs_modulus(Real E) {
        E_ = E;
        if (is_2d_) {
            compute_D_matrix_2D();
        } else {
            compute_D_matrix_3D();
        }
    }

    void set_poissons_ratio(Real nu) {
        nu_ = nu;
        if (is_2d_) {
            compute_D_matrix_2D();
        } else {
            compute_D_matrix_3D();
        }
    }

    void set_plane_type(PlaneType type) {
        plane_type_ = type;
        if (is_2d_) {
            compute_D_matrix_2D();
        }
    }

    Real youngs_modulus() const { return E_; }
    Real poissons_ratio() const { return nu_; }
    PlaneType plane_type() const { return plane_type_; }
    bool is_2d() const { return is_2d_; }

private:
    Real E_;                ///< 杨氏模量
    Real nu_;               ///< 泊松比
    bool is_2d_;            ///< 2D 模式开关
    PlaneType plane_type_;  ///< 平面类型 (仅 2D)

    DenseMatrix D_;         ///< 本构矩阵 (3x3 for 2D, 6x6 for 3D)

    /**
     * 计算 2D 本构矩阵 D
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
    void compute_D_matrix_2D();
    
    /**
     * 计算 3D 本构矩阵 D
     * 
     * D = E/((1+ν)(1-2ν)) * [[1-ν,  ν,    ν,    0,         0,         0        ],
     *                         [ν,    1-ν,  ν,    0,         0,         0        ],
     *                         [ν,    ν,    1-ν,  0,         0,         0        ],
     *                         [0,    0,    0,    (1-2ν)/2,  0,         0        ],
     *                         [0,    0,    0,    0,         (1-2ν)/2,  0        ],
     *                         [0,    0,    0,    0,         0,         (1-2ν)/2 ]]
     */
    void compute_D_matrix_3D();
};

} // namespace physics
} // namespace fem
