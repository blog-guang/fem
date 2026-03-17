/**
 * dof_handler.h - 自由度管理器
 * 
 * 负责:
 * - DOF 编号分配
 * - 稀疏模式构建
 * - local DOF ↔ global DOF 映射
 */

#pragma once

#include "core/types.h"
#include "mesh/model.h"
#include "math/sparse_matrix.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

namespace fem {

// 场类型定义 (使用 core/types.h 中的 FieldType)
// DISPLACEMENT = Vector3D, PRESSURE = Scalar, TEMPERATURE = Scalar

// 有限元空间描述
struct FiniteElementSpace {
    FieldType field_type;     // 场类型
    int order;                // 单元阶次 (1=线性, 2=二次)
    int components;           // 分量数量 (位移: 3, 标量: 1)
    
    FiniteElementSpace(FieldType ft = FieldType::DISPLACEMENT, 
                       int o = 1, int c = 3)
        : field_type(ft), order(o), components(c) {}
};

/**
 * DOF 索引信息
 */
struct DofIndex {
    Index element_id;    // 所在单元 ID
    int local_dof;       // 单元内局部 DOF 编号
    Index global_dof;    // 全局 DOF 编号
    
    bool operator<(const DofIndex& other) const {
        return global_dof < other.global_dof;
    }
};

/**
 * DofHandler - 自由度管理器
 * 
 * 使用方法:
 * ```cpp
 * DofHandler dof_handler(model);
 * dof_handler.distribute_dofs(fe_space);
 * 
 * auto pattern = dof_handler.make_sparsity_pattern();
 * auto K = SparseMatrix(pattern);
 * ```
 */
class DofHandler {
public:
    /**
     * 构造函数
     * @param model 有限元模型
     */
    explicit DofHandler(const Model& model);
    
    // ═══ DOF 编号 ═══
    
    /**
     * 分配全局 DOF 编号
     * @param fe_space 有限元空间描述
     */
    void distribute_dofs(const FiniteElementSpace& fe_space);
    
    /**
     * 获取总 DOF 数量
     */
    Index n_dofs() const { return n_dofs_; }
    
    /**
     * 获取每节点 DOF 数
     */
    Index dofs_per_node() const { return dofs_per_node_; }
    
    // ═══ DOF 映射 ═══
    
    /**
     * 获取节点的所有 DOF 全局编号
     * @param node_id 节点 ID
     * @return DOF 全局编号列表
     */
    std::vector<Index> node_dofs(Index node_id) const;
    
    /**
     * 获取单元的所有 DOF 全局编号
     * @param elem_id 单元 ID
     * @return DOF 全局编号列表
     */
    std::vector<Index> element_dofs(Index elem_id) const;
    
    /**
     * 通过节点 ID 和分量获取 DOF 编号
     */
    Index dof_id(Index node_id, int component) const;
    
    // ═══ 稀疏模式 ═══
    
    /**
     * 构建稀疏模式（两遍扫描）
     * @return 稀疏矩阵模式
     */
    SparseMatrixPattern make_sparsity_pattern() const;
    
    /**
     * 估计非零元数量
     */
    size_t estimate_nnz() const;
    
    // ═══ 边界条件支持 ═══
    
    /**
     * 标记 Dirichlet 边界条件 DOF
     * @param node_id 节点 ID
     * @param component 分量 (0=x, 1=y, 2=z)
     * @param value 边界值
     */
    void add_dirichlet_bc(Index node_id, int component, Real value);
    
    /**
     * 获取所有 Dirichlet DOF
     */
    const std::map<Index, Real>& dirichlet_dofs() const { return dirichlet_bcs_; }
    
    /**
     * 检查是否是 Dirichlet DOF
     */
    bool is_dirichlet_dof(Index global_dof) const {
        return dirichlet_bcs_.find(global_dof) != dirichlet_bcs_.end();
    }
    
    // ═══ 信息输出 ═══
    void print_info() const;
    
private:
    const Model& model_;                    ///< 模型引用
    FiniteElementSpace fe_space_;           ///< 有限元空间
    
    Index n_dofs_ = 0;                     ///< 总 DOF 数
    Index dofs_per_node_ = 1;              ///< 每节点 DOF 数
    int mesh_id_ = 0;                      ///< 使用的网格 ID
    
    // DOF 映射表: node_id * dofs_per_node + component → global_dof
    std::vector<Index> node_to_dof_;      
    
    // 单元 DOF 缓存
    mutable std::map<Index, std::vector<Index>> element_dof_cache_;
    
    // Dirichlet 边界条件
    std::map<Index, Real> dirichlet_bcs_;
    
    // ═══ 内部方法 ═══
    
    /**
     * 计算每节点 DOF 数
     */
    void compute_dofs_per_node();
    
    /**
     * 构建节点到 DOF 的映射
     */
    void build_node_dof_map();
    
    /**
     * 构建单元 DOF 缓存
     */
    void build_element_dof_cache() const;
};

/**
 * 稀疏矩阵模式构建器
 * 
 * 用于预分配稀疏矩阵内存
 */
class SparsityPatternBuilder {
public:
    SparsityPatternBuilder(Index n_dofs);
    
    /**
     * 添加非零元位置
     */
    void add(Index row, Index col);
    
    /**
     * 添加多个 DOF 之间的连接（单元矩阵）
     */
    void add_dofs(const std::vector<Index>& dofs);
    
    /**
     * 构建最终模式
     */
    SparseMatrixPattern build() const;
    
    /**
     * 获取估计的非零元数量
     */
    size_t nnz() const { return entries_.size(); }
    
private:
    Index n_dofs_;
    std::vector<std::pair<Index, Index>> entries_;
};

}  // namespace fem
