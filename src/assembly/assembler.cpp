#include "assembly/assembler.h"
#include "core/logger.h"

namespace fem {

Assembler::Assembler(const Model& model, Index dofs_per_node)
    : model_(model), dofs_per_node_(dofs_per_node) {
    
    if (dofs_per_node_ == 0) {
        throw std::invalid_argument("dofs_per_node must be > 0");
    }

    // 计算总自由度数 (假设所有网格节点数相同)
    n_dofs_ = 0;
    for (std::size_t i = 0; i < model_.num_meshes(); ++i) {
        const Mesh& mesh = model_.mesh(i);
        n_dofs_ = std::max(n_dofs_, mesh.num_nodes() * dofs_per_node_);
    }

    if (n_dofs_ == 0) {
        throw std::runtime_error("Model has no nodes");
    }

    // 初始化 COO 矩阵和向量
    K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);
    F_ = Vector(n_dofs_, 0.0);

    FEM_INFO("Assembler initialized: " + std::to_string(n_dofs_) + " DOFs (" +
             std::to_string(dofs_per_node_) + " per node)");
}

void Assembler::assemble(ElementMatrixFunc elem_func) {
    clear();

    // 遍历所有网格
    for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
        const Mesh& mesh = model_.mesh(mesh_id);

        // 遍历所有单元
        for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            const Element& elem = mesh.element(elem_id);
            
            Index n_nodes = elem.nodes().size();
            Index n_dofs_elem = n_nodes * dofs_per_node_;

            // 单元矩阵和向量
            DenseMatrix Ke(n_dofs_elem, n_dofs_elem, 0.0);
            Vector Fe(n_dofs_elem, 0.0);

            // 调用用户提供的单元计算函数
            elem_func(elem_id, mesh, Ke, Fe);

            // 装配到全局系统
            for (Index i = 0; i < n_nodes; ++i) {
                Index node_i = elem.nodes()[i];

                for (Index di = 0; di < dofs_per_node_; ++di) {
                    Index gi = node_i * dofs_per_node_ + di;  // 全局DOF索引

                    // 装配向量
                    F_[gi] += Fe[i * dofs_per_node_ + di];

                    // 装配矩阵
                    for (Index j = 0; j < n_nodes; ++j) {
                        Index node_j = elem.nodes()[j];

                        for (Index dj = 0; dj < dofs_per_node_; ++dj) {
                            Index gj = node_j * dofs_per_node_ + dj;
                            
                            K_coo_.add(gi, gj, Ke(i * dofs_per_node_ + di, 
                                                   j * dofs_per_node_ + dj));
                        }
                    }
                }
            }
        }
    }

    assembled_ = true;
    FEM_INFO("Assembly completed: " + std::to_string(K_coo_.nnz()) + " non-zeros");
}

void Assembler::assemble_matrix(ElementMatrixFunc elem_func) {
    K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);

    for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
        const Mesh& mesh = model_.mesh(mesh_id);

        for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            const Element& elem = mesh.element(elem_id);
            
            Index n_nodes = elem.nodes().size();
            Index n_dofs_elem = n_nodes * dofs_per_node_;

            DenseMatrix Ke(n_dofs_elem, n_dofs_elem, 0.0);
            Vector Fe(n_dofs_elem, 0.0);  // 占位

            elem_func(elem_id, mesh, Ke, Fe);

            // 装配矩阵
            for (Index i = 0; i < n_nodes; ++i) {
                Index node_i = elem.nodes()[i];

                for (Index di = 0; di < dofs_per_node_; ++di) {
                    Index gi = node_i * dofs_per_node_ + di;

                    for (Index j = 0; j < n_nodes; ++j) {
                        Index node_j = elem.nodes()[j];

                        for (Index dj = 0; dj < dofs_per_node_; ++dj) {
                            Index gj = node_j * dofs_per_node_ + dj;
                            
                            K_coo_.add(gi, gj, Ke(i * dofs_per_node_ + di, 
                                                   j * dofs_per_node_ + dj));
                        }
                    }
                }
            }
        }
    }

    assembled_ = true;
}

void Assembler::assemble_vector(ElementMatrixFunc elem_func) {
    F_ = Vector(n_dofs_, 0.0);

    for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
        const Mesh& mesh = model_.mesh(mesh_id);

        for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            const Element& elem = mesh.element(elem_id);
            
            Index n_nodes = elem.nodes().size();
            Index n_dofs_elem = n_nodes * dofs_per_node_;

            DenseMatrix Ke(n_dofs_elem, n_dofs_elem, 0.0);  // 占位
            Vector Fe(n_dofs_elem, 0.0);

            elem_func(elem_id, mesh, Ke, Fe);

            // 装配向量
            for (Index i = 0; i < n_nodes; ++i) {
                Index node_i = elem.nodes()[i];

                for (Index di = 0; di < dofs_per_node_; ++di) {
                    Index gi = node_i * dofs_per_node_ + di;
                    F_[gi] += Fe[i * dofs_per_node_ + di];
                }
            }
        }
    }
}

void Assembler::apply_dirichlet(const std::vector<DirichletBC>& bcs) {
    if (!assembled_) {
        throw std::runtime_error("Must call assemble() before apply_dirichlet()");
    }

    // 转换为 CSR 格式以便修改
    SparseMatrixCSR K_csr = coo_to_csr(K_coo_);

    // 收集所有需要约束的 DOF 及其值
    std::vector<Real> bc_values(n_dofs_, 0.0);
    std::vector<bool> is_bc_dof(n_dofs_, false);

    for (const auto& bc : bcs) {
        for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
            const Mesh& mesh = model_.mesh(mesh_id);

            if (!mesh.has_boundary(bc.boundary_name)) {
                continue;
            }

            const auto& boundary_nodes = mesh.boundary(bc.boundary_name);

            for (Index node_id : boundary_nodes) {
                Index dof_id = node_id * dofs_per_node_ + bc.dof;
                if (dof_id < n_dofs_) {
                    is_bc_dof[dof_id] = true;
                    bc_values[dof_id] = bc.value;
                }
            }
        }
    }

    // 统一应用边界条件 (简化方法: 主对角线法)
    // K(ii) = 1, K(i,j≠i) = 0, F(i) = value
    for (Index i = 0; i < n_dofs_; ++i) {
        if (is_bc_dof[i]) {
            Index row_start = K_csr.row_ptr()[i];
            Index row_end = K_csr.row_ptr()[i + 1];

            for (Index k = row_start; k < row_end; ++k) {
                Index j = K_csr.col_indices()[k];
                if (j == i) {
                    K_csr.values()[k] = 1.0;
                } else {
                    K_csr.values()[k] = 0.0;
                }
            }

            F_[i] = bc_values[i];
        }
    }

    // 转换回 COO 格式
    K_coo_ = csr_to_coo(K_csr);

    FEM_INFO("Applied " + std::to_string(bcs.size()) + " Dirichlet BCs");
}

SparseMatrixCSR Assembler::matrix() const {
    if (!assembled_) {
        throw std::runtime_error("Must call assemble() before matrix()");
    }
    return coo_to_csr(K_coo_);
}

void Assembler::clear() {
    K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);
    F_ = Vector(n_dofs_, 0.0);
    assembled_ = false;
}

} // namespace fem
