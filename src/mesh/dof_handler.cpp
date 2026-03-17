/**
 * dof_handler.cpp - 自由度管理器实现 (Refactored)
 */

#include "mesh/dof_handler.h"
#include "core/logger.h"
#include <algorithm>

namespace fem {

// ═══════════════════════════════════════════════════════════
// DofHandler 实现
// ═══════════════════════════════════════════════════════════

DofHandler::DofHandler(const Model& model, int mesh_id)
    : model_(model), mesh_id_(mesh_id)
{
    if (mesh_id_ < 0 || mesh_id_ >= static_cast<int>(model_.num_meshes())) {
        throw std::out_of_range("Invalid mesh_id: " + std::to_string(mesh_id_));
    }
}

void DofHandler::distribute_dofs(const FiniteElementSpace& fe_space) {
    fe_space_ = fe_space;
    dofs_per_node_ = fe_space.components();
    
    const Mesh& mesh = model_.mesh(mesh_id_);
    n_dofs_ = mesh.num_nodes() * dofs_per_node_;
    
    FEM_INFO("DofHandler: " + std::to_string(n_dofs_) + " DOFs, " +
             std::to_string(dofs_per_node_) + " DOFs/node, " +
             std::to_string(mesh.num_nodes()) + " nodes");
}

size_t DofHandler::element_dofs(Index elem_id, std::vector<Index>& dofs) const {
    const Mesh& mesh = model_.mesh(mesh_id_);
    
    if (elem_id >= mesh.num_elements()) {
        throw std::out_of_range("Invalid elem_id: " + std::to_string(elem_id));
    }
    
    const auto& elem = mesh.element(elem_id);
    const auto& node_ids = elem.nodes();
    
    size_t num_elem_dofs = node_ids.size() * dofs_per_node_;
    dofs.resize(num_elem_dofs);
    
    size_t idx = 0;
    for (Index node_id : node_ids) {
        for (int comp = 0; comp < dofs_per_node_; comp++) {
            dofs[idx++] = node_dof(node_id, comp);
        }
    }
    
    return num_elem_dofs;
}

std::vector<Index> DofHandler::element_dofs(Index elem_id) const {
    std::vector<Index> dofs;
    element_dofs(elem_id, dofs);
    return dofs;
}

SparseMatrixPattern DofHandler::make_sparsity_pattern() const {
    SparsityPatternBuilder builder(n_dofs_);
    
    const Mesh& mesh = model_.mesh(mesh_id_);
    std::vector<Index> elem_dofs;
    
    // 遍历所有单元，添加 DOF 连接
    for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
        element_dofs(elem_id, elem_dofs);
        builder.add_element_dofs(elem_dofs);
    }
    
    return builder.build();
}

size_t DofHandler::estimate_nnz() const {
    const Mesh& mesh = model_.mesh(mesh_id_);
    
    // 保守估计：每个单元的 DOF 对 
    // (实际会因为节点共享而减少)
    size_t avg_nodes_per_elem = 4;  // 假设 Tet4/Quad4
    size_t dofs_per_elem = avg_nodes_per_elem * dofs_per_node_;
    size_t nnz_per_elem = dofs_per_elem * dofs_per_elem;
    
    return mesh.num_elements() * nnz_per_elem;
}

void DofHandler::print_info() const {
    FEM_INFO("=== DofHandler Information ===");
    FEM_INFO("  Total DOFs: " + std::to_string(n_dofs_));
    FEM_INFO("  DOFs per node: " + std::to_string(dofs_per_node_));
    FEM_INFO("  Mesh ID: " + std::to_string(mesh_id_));
    FEM_INFO("  Field type: " + std::to_string(static_cast<int>(fe_space_.field_type)));
    FEM_INFO("  Order: " + std::to_string(fe_space_.order));
}

// ═══════════════════════════════════════════════════════════
// SparsityPatternBuilder 实现 (优化版)
// ═══════════════════════════════════════════════════════════

SparsityPatternBuilder::SparsityPatternBuilder(Index n_dofs)
    : n_dofs_(n_dofs)
{
    row_entries_.resize(n_dofs_);
}

void SparsityPatternBuilder::add_element_dofs(const std::vector<Index>& dofs) {
    // 单元内所有 DOF 两两连接
    for (Index i : dofs) {
        if (i < n_dofs_) {
            for (Index j : dofs) {
                if (j < n_dofs_) {
                    row_entries_[i].insert(j);  // set 自动去重和排序
                }
            }
        }
    }
    
    nnz_estimate_ += dofs.size() * dofs.size();
}

SparseMatrixPattern SparsityPatternBuilder::build() const {
    // 1. 计算 row_ptr
    std::vector<size_t> row_ptr(n_dofs_ + 1, 0);
    
    for (Index i = 0; i < n_dofs_; i++) {
        row_ptr[i + 1] = row_ptr[i] + row_entries_[i].size();
    }
    
    size_t nnz = row_ptr[n_dofs_];
    
    // 2. 构建 col_indices
    std::vector<Index> col_indices;
    col_indices.reserve(nnz);
    
    for (Index i = 0; i < n_dofs_; i++) {
        // set 已经排序，直接复制
        col_indices.insert(col_indices.end(), 
                          row_entries_[i].begin(), 
                          row_entries_[i].end());
    }
    
    return SparseMatrixPattern(n_dofs_, n_dofs_, 
                               std::move(row_ptr), 
                               std::move(col_indices));
}

}  // namespace fem
