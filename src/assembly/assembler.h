/**
 * assembler.h - 通用有限元装配器
 * 
 * 支持:
 * - 刚度矩阵装配 (K)
 * - 载荷向量装配 (F)
 * - Dirichlet 边界条件
 * - 多自由度场 (标量/矢量)
 */

#pragma once

#include "core/types.h"
#include "mesh/model.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "math/sparse_matrix.h"

#include <functional>
#include <vector>

namespace fem {

/**
 * 单元矩阵/向量计算函数类型
 * 
 * @param elem_id 单元ID
 * @param mesh 网格引用
 * @param Ke 输出：单元矩阵 (nnode*dof x nnode*dof)
 * @param Fe 输出：单元向量 (nnode*dof)
 */
using ElementMatrixFunc = std::function<void(
    Index elem_id, 
    const Mesh& mesh,
    DenseMatrix& Ke,
    Vector& Fe
)>;

/**
 * 边界条件数据
 */
struct DirichletBC {
    std::string boundary_name;  ///< 边界名称 (e.g., "left", "right")
    Index dof;                  ///< 自由度索引 (0=x, 1=y, etc.)
    Real value;                 ///< 边界值
};

/**
 * Assembler - 有限元装配器
 * 
 * 使用方法:
 * ```cpp
 * Assembler asm(model, dofs_per_node);
 * asm.assemble(element_matrix_func);
 * asm.apply_dirichlet(bcs);
 * auto K = asm.matrix();
 * auto F = asm.rhs();
 * ```
 */
class Assembler {
public:
    /**
     * 构造函数
     * @param model 模型引用
     * @param dofs_per_node 每节点自由度数 (1=标量, 2=2D矢量, 3=3D矢量)
     */
    Assembler(const Model& model, Index dofs_per_node = 1);

    /**
     * 装配全局矩阵和向量
     * @param elem_func 单元矩阵/向量计算函数
     */
    void assemble(ElementMatrixFunc elem_func);

    /**
     * 仅装配全局矩阵
     * @param elem_func 单元矩阵计算函数 (Fe 可忽略)
     */
    void assemble_matrix(ElementMatrixFunc elem_func);

    /**
     * 仅装配全局向量
     * @param elem_func 单元向量计算函数 (Ke 可忽略)
     */
    void assemble_vector(ElementMatrixFunc elem_func);

    /**
     * 应用 Dirichlet 边界条件
     * 
     * 方法: 主对角线置1法
     * - K(ii) = 1.0
     * - K(i,j) = 0.0 (j≠i)
     * - F(i) = bc_value
     * 
     * @param bcs 边界条件列表
     */
    void apply_dirichlet(const std::vector<DirichletBC>& bcs);

    /**
     * 获取装配后的刚度矩阵 (CSR格式)
     */
    SparseMatrixCSR matrix() const;

    /**
     * 获取装配后的载荷向量
     */
    const Vector& rhs() const { return F_; }

    /**
     * 清空矩阵和向量
     */
    void clear();

    /**
     * 获取总自由度数
     */
    Index num_dofs() const { return n_dofs_; }

    /**
     * 获取每节点自由度数
     */
    Index dofs_per_node() const { return dofs_per_node_; }

private:
    const Model& model_;          ///< 模型引用
    Index dofs_per_node_;         ///< 每节点自由度数
    Index n_dofs_;                ///< 总自由度数

    SparseMatrixCOO K_coo_;       ///< 刚度矩阵 (COO格式，装配阶段)
    Vector F_;                    ///< 载荷向量

    bool assembled_{false};       ///< 是否已装配
};

} // namespace fem
