/**
 * physics_base.h - 物理模块基类
 * 
 * 提供统一的单元装配框架：
 * - 形函数系统集成
 * - 高斯积分循环
 * - B 矩阵构造（应变-位移关系）
 * - 支持所有单元类型（2D/3D）
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "math/dense_matrix.h"
#include "math/vector.h"
#include "shape/shape_function_factory.h"

namespace fem {
namespace physics {

/**
 * PhysicsBase - 物理模块基类
 * 
 * 提供通用的单元装配工具：
 * - 获取节点坐标
 * - 高斯积分循环
 * - 形函数导数计算
 * - B 矩阵构造（弹性力学）
 * 
 * 子类只需实现具体的物理方程。
 */
class PhysicsBase {
public:
    virtual ~PhysicsBase() = default;
    
    /**
     * 计算单元刚度矩阵和载荷向量（纯虚函数）
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @param Ke 输出：单元刚度矩阵
     * @param Fe 输出：单元载荷向量
     */
    virtual void compute_element(Index elem_id, const Mesh& mesh,
                                DenseMatrix& Ke, Vector& Fe) const = 0;

protected:
    // ═══ 通用工具函数 ═══
    
    /**
     * 获取单元节点坐标
     * 
     * @param elem_id 单元ID
     * @param mesh 网格引用
     * @return 节点坐标向量
     */
    std::vector<Vec3> get_element_coords(Index elem_id, const Mesh& mesh) const;
    
    /**
     * 检查单元维度是否匹配
     * 
     * @param elem 单元引用
     * @param expected_dim 期望维度 (2 or 3)
     * @return true 如果匹配
     */
    bool check_dimension(const Element& elem, int expected_dim) const;
    
    /**
     * 计算单元维度和自由度数
     * 
     * @param elem 单元引用
     * @param dofs_per_node 每节点自由度数
     * @return [n_nodes, n_dofs, dim]
     */
    std::tuple<Index, Index, int> get_element_info(
        const Element& elem, int dofs_per_node) const;
    
    // ═══ 高斯积分工具 ═══
    
    /**
     * 为单元类型选择高斯积分阶数
     * 
     * @param elem_type 单元类型
     * @return 积分阶数（默认2阶）
     */
    int get_gauss_order(ElementType elem_type) const;
    
    // ═══ B 矩阵构造（弹性力学） ═══
    
    /**
     * 构造 2D 弹性力学 B 矩阵
     * 
     * B = [[dN1/dx, 0,      dN2/dx, 0,      ...],
     *      [0,      dN1/dy, 0,      dN2/dy, ...],
     *      [dN1/dy, dN1/dx, dN2/dy, dN2/dx, ...]]
     * 
     * @param dN_dxyz 形函数物理坐标导数 (n_nodes × 2)
     * @return B 矩阵 (3 × n_dofs)
     */
    DenseMatrix build_B_matrix_2D(const DenseMatrix& dN_dxyz) const;
    
    /**
     * 构造 3D 弹性力学 B 矩阵
     * 
     * B = [[dN1/dx, 0,      0,      ...],
     *      [0,      dN1/dy, 0,      ...],
     *      [0,      0,      dN1/dz, ...],
     *      [0,      dN1/dz, dN1/dy, ...],  // γ_yz
     *      [dN1/dz, 0,      dN1/dx, ...],  // γ_xz
     *      [dN1/dy, dN1/dx, 0,      ...]]  // γ_xy
     * 
     * @param dN_dxyz 形函数物理坐标导数 (n_nodes × 3)
     * @return B 矩阵 (6 × n_dofs)
     */
    DenseMatrix build_B_matrix_3D(const DenseMatrix& dN_dxyz) const;
    
    // ═══ 雅可比计算 ═══
    
    /**
     * 计算雅可比行列式
     * 
     * @param J 雅可比矩阵
     * @return det(J)
     */
    Real compute_jacobian_determinant(const DenseMatrix& J) const;
};

} // namespace physics
} // namespace fem
