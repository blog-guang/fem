/**
 * dof_handler.cpp - 自由度管理器实现
 */

#include "mesh/dof_handler.h"
#include "core/logger.h"
#include <algorithm>
#include <cassert>

namespace fem {

// ═══════════════════════════════════════════════════════════
// DofHandler 实现
// ═══════════════════════════════════════════════════════════

DofHandler::DofHandler(const Model& model)
    : model_(model)
{
}

void DofHandler::distribute_dofs(const FiniteElementSpace& fe_space) {
    fe_space_ = fe_space;
    
    // 计算每节点 DOF 数
    compute_dofs_per_node();
    
    const Mesh& mesh = model_.mesh(mesh_id_);
    Index num_nodes = mesh.num_nodes();
    
    // 分配 DOF 编号: node_id * dofs_per_node + component
    n_dofs_ = num_nodes * dofs_per_node_;
    node_to_dof_.resize(num_nodes * dofs_per_node_);
    
    for (Index node_id = 0; node_id < num_nodes; node_id++) {
        for (int comp = 0; comp < dofs_per_node_; comp++) {
            Index global_dof = node_id * dofs_per_node_ + comp;
            node_to_dof_[node_id * dofs_per_node_ + comp] = global_dof;
        }
    }
    
    // 构建单元 DOF 缓存
    build_element_dof_cache();
    
    FEM_INFO("DofHandler: " + std::to_string(n_dofs_) + " DOFs, " +
             std::to_string(dofs_per_node_) + " DOFs/node");
}

void DofHandler::compute_dofs_per_node() {
    switch (fe_space_.field_type) {
        case FieldType::Vector3D:
            dofs_per_node_ = 3;
            break;
        case FieldType::Vector2D:
            dofs_per_node_ = 2;
            break;
        case FieldType::Scalar:
        default:
            dofs_per_node_ = 1;
            break;
    }
}

void DofHandler::build_node_dof_map() {
    const Mesh& mesh = model_.mesh(mesh_id_);
    Index num_nodes = mesh.num_nodes();
    
    node_to_dof_.resize(num_nodes * dofs_per_node_);
    
    for (Index node_id = 0; node_id < num_nodes; node_id++) {
        for (int comp = 0; comp < dofs_per_node_; comp++) {
            node_to_dof_[node_id * dofs_per_node_ + comp] = 
                node_id * dofs_per_node_ + comp;
        }
    }
}

void DofHandler::build_element_dof_cache() const {
    const Mesh& mesh = model_.mesh(mesh_id_);
    element_dof_cache_.clear();
    
    for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
        const auto& elem = mesh.element(elem_id);
        const auto& node_ids = elem.nodes();
        
        std::vector<Index> dofs;
        dofs.reserve(node_ids.size() * dofs_per_node_);
        
        for (Index node_id : node_ids) {
            for (int comp = 0; comp < dofs_per_node_; comp++) {
                dofs.push_back(node_id * dofs_per_node_ + comp);
            }
        }
        
        element_dof_cache_[elem_id] = std::move(dofs);
    }
}

std::vector<Index> DofHandler::node_dofs(Index node_id) const {
    std::vector<Index> dofs;
    dofs.reserve(dofs_per_node_);
    
    for (int comp = 0; comp < dofs_per_node_; comp++) {
        dofs.push_back(node_id * dofs_per_node_ + comp);
    }
    
    return dofs;
}

std::vector<Index> DofHandler::element_dofs(Index elem_id) const {
    auto it = element_dof_cache_.find(elem_id);
    if (it != element_dof_cache_.end()) {
        return it->second;
    }
    
    // 如果缓存中没有，重新计算
    const Mesh& mesh = model_.mesh(mesh_id_);
    const auto& elem = mesh.element(elem_id);
    const auto& node_ids = elem.nodes();
    
    std::vector<Index> dofs;
    dofs.reserve(node_ids.size() * dofs_per_node_);
    
    for (Index node_id : node_ids) {
        for (int comp = 0; comp < dofs_per_node_; comp++) {
            dofs.push_back(node_id * dofs_per_node_ + comp);
        }
    }
    
    return dofs;
}

Index DofHandler::dof_id(Index node_id, int component) const {
    assert(node_id * dofs_per_node_ + component < node_to_dof_.size());
    return node_to_dof_[node_id * dofs_per_node_ + component];
}

SparseMatrixPattern DofHandler::make_sparsity_pattern() const {
    SparsityPatternBuilder builder(n_dofs_);
    
    const Mesh& mesh = model_.mesh(mesh_id_);
    
    // 遍历所有单元，添加 DOF 连接
    for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
        auto dofs = element_dofs(elem_id);
        builder.add_dofs(dofs);
    }
    
    return builder.build();
}

size_t DofHandler::estimate_nnz() const {
    const Mesh& mesh = model_.mesh(mesh_id_);
    
    // 估计：每个单元的 DOF 对 connections
    // 对于 Lagrange 单元：(dofs_per_elem)^2
    size_t nnz_per_elem = dofs_per_node_ * dofs_per_node_;
    
    // 实际非零元会比这个少（因为共享节点），但保守估计
    return mesh.num_elements() * nnz_per_elem;
}

void DofHandler::add_dirichlet_bc(Index node_id, int component, Real value) {
    Index global_dof = dof_id(node_id, component);
    dirichlet_bcs_[global_dof] = value;
}

void DofHandler::print_info() const {
    FEM_INFO("=== DofHandler Information ===");
    FEM_INFO("Total DOFs: " + std::to_string(n_dofs_));
    FEM_INFO("DOFs per node: " + std::to_string(dofs_per_node_));
    FEM_INFO("Mesh ID: " + std::to_string(mesh_id_));
    FEM_INFO("Dirichlet BCs: " + std::to_string(dirichlet_bcs_.size()));
}

// ═══════════════════════════════════════════════════════════
// SparsityPatternBuilder 实现
// ═══════════════════════════════════════════════════════════

SparsityPatternBuilder::SparsityPatternBuilder(Index n_dofs)
    : n_dofs_(n_dofs)
{
    entries_.reserve(n_dofs * 10);  // 预分配
}

void SparsityPatternBuilder::add(Index row, Index col) {
    if (row < n_dofs_ && col < n_dofs_) {
        entries_.push_back({row, col});
    }
}

void SparsityPatternBuilder::add_dofs(const std::vector<Index>& dofs) {
    // 单元内所有 DOF 两两连接
    for (size_t i = 0; i < dofs.size(); i++) {
        for (size_t j = 0; j < dofs.size(); j++) {
            add(dofs[i], dofs[j]);
        }
    }
}

SparseMatrixPattern SparsityPatternBuilder::build() const {
    // 1. 排序并去重
    std::vector<std::pair<Index, Index>> unique_entries = entries_;
    std::sort(unique_entries.begin(), unique_entries.end());
    unique_entries.erase(
        std::unique(unique_entries.begin(), unique_entries.end()),
        unique_entries.end()
    );
    
    // 2. 统计每行非零元数量
    std::vector<size_t> nnz_per_row(n_dofs_, 0);
    for (const auto& [row, col] : unique_entries) {
        nnz_per_row[row]++;
    }
    
    // 3. 构建行指针
    std::vector<size_t> row_ptr(n_dofs_ + 1, 0);
    for (Index i = 0; i < n_dofs_; i++) {
        row_ptr[i + 1] = row_ptr[i] + nnz_per_row[i];
    }
    
    // 4. 构建列索引
    std::vector<Index> col_idx(unique_entries.size());
    std::vector<size_t> current_pos(n_dofs_, 0);
    
    for (const auto& [row, col] : unique_entries) {
        size_t idx = row_ptr[row] + current_pos[row]++;
        col_idx[idx] = col;
    }
    
    return SparseMatrixPattern(n_dofs_, n_dofs_, 
                               std::move(row_ptr), 
                               std::move(col_idx));
}

}  // namespace fem
