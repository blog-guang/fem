/**
 * dof_handler.h - 自由度管理器 (Refactored)
 * 
 * 职责:
 * - DOF 编号分配与映射
 * - 稀疏模式构建
 * - 提供高效的 DOF 查询接口
 */

#pragma once

#include "core/types.h"
#include "mesh/model.h"
#include "math/sparse_matrix.h"
#include <vector>
#include <set>
#include <memory>

namespace fem {

/**
 * 有限元空间描述
 */
struct FiniteElementSpace {
    FieldType field_type;     ///< 场类型
    int order;                ///< 单元阶次 (1=线性, 2=二次)
    
    FiniteElementSpace(FieldType ft = FieldType::Vector3D, int o = 1)
        : field_type(ft), order(o) {}
    
    /// 获取每节点分量数
    int components() const {
        switch (field_type) {
            case FieldType::Vector3D: return 3;
            case FieldType::Vector2D: return 2;
            case FieldType::Scalar:
            default: return 1;
        }
    }
};

/**
 * DofHandler - 自由度管理器
 * 
 * 设计原则:
 * - 轻量级：不存储冗余数据
 * - 高性能：内联计算 DOF 编号
 * - 线程安全：无可变状态
 * 
 * 使用方法:
 * ```cpp
 * DofHandler dof_handler(model);
 * dof_handler.distribute_dofs(fe_space);
 * 
 * // 查询 DOF
 * Index dof = dof_handler.node_dof(node_id, component);
 * 
 * // 构建稀疏模式
 * auto pattern = dof_handler.make_sparsity_pattern();
 * ```
 */
class DofHandler {
public:
    /**
     * 构造函数
     * @param model 有限元模型
     * @param mesh_id 使用的网格 ID (默认 0)
     */
    explicit DofHandler(const Model& model, int mesh_id = 0);
    
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
    int dofs_per_node() const { return dofs_per_node_; }
    
    // ═══ DOF 映射 (高性能内联接口) ═══
    
    /**
     * 获取节点某分量的全局 DOF 编号
     * @param node_id 节点 ID
     * @param component 分量索引 (0-based)
     * @return 全局 DOF 编号
     */
    Index node_dof(Index node_id, int component) const {
        return node_id * dofs_per_node_ + component;
    }
    
    /**
     * 获取单元的所有 DOF (预分配 buffer 版本，避免内存分配)
     * @param elem_id 单元 ID
     * @param dofs 输出缓冲区 (应预分配足够空间)
     * @return 实际 DOF 数量
     */
    size_t element_dofs(Index elem_id, std::vector<Index>& dofs) const;
    
    /**
     * 获取单元的所有 DOF (便捷版本)
     * @param elem_id 单元 ID
     * @return DOF 索引向量
     */
    std::vector<Index> element_dofs(Index elem_id) const;
    
    // ═══ 稀疏模式 ═══
    
    /**
     * 构建稀疏模式
     * @return 稀疏矩阵模式
     */
    SparseMatrixPattern make_sparsity_pattern() const;
    
    /**
     * 估计非零元数量 (用于预分配)
     */
    size_t estimate_nnz() const;
    
    // ═══ 信息查询 ═══
    
    const Model& model() const { return model_; }
    int mesh_id() const { return mesh_id_; }
    const FiniteElementSpace& fe_space() const { return fe_space_; }
    
    void print_info() const;
    
private:
    const Model& model_;                    ///< 模型引用
    int mesh_id_;                           ///< 使用的网格 ID
    FiniteElementSpace fe_space_;           ///< 有限元空间
    
    Index n_dofs_ = 0;                     ///< 总 DOF 数
    int dofs_per_node_ = 1;                ///< 每节点 DOF 数
};

/**
 * 稀疏矩阵模式构建器 (优化版)
 * 
 * 使用 set 自动去重，避免排序开销
 */
class SparsityPatternBuilder {
public:
    explicit SparsityPatternBuilder(Index n_dofs);
    
    /**
     * 添加单元的所有 DOF 连接
     * @param dofs 单元 DOF 索引列表
     */
    void add_element_dofs(const std::vector<Index>& dofs);
    
    /**
     * 构建最终模式 (CSR 格式)
     */
    SparseMatrixPattern build() const;
    
    /**
     * 获取当前非零元数量 (估计值，可能有重复)
     */
    size_t nnz_estimate() const { return nnz_estimate_; }
    
private:
    Index n_dofs_;
    
    // 使用 set 按行存储列索引，自动去重和排序
    std::vector<std::set<Index>> row_entries_;
    
    size_t nnz_estimate_ = 0;  ///< 非零元估计值
};

}  // namespace fem
