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
    is_dirichlet_dof_.resize(n_dofs_, false);

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
    
    // 保存原始刚度矩阵用于反力计算
    if (!has_original_) {
        K_original_ = K_csr;
        has_original_ = true;
    }

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
                    is_dirichlet_dof_[dof_id] = true;  // 标记为 Dirichlet DOF
                }
            }
        }
    }

    // 完全消去法 (保持系统一致性)
    // 1. 先修正右端项: F(j) -= K(j,i) * bc_value (对所有自由DOF j)
    for (Index i = 0; i < n_dofs_; ++i) {
        if (is_bc_dof[i]) {
            Real bc_val = bc_values[i];
            
            // 遍历所有行，找到列i的元素
            for (Index row = 0; row < n_dofs_; ++row) {
                if (is_bc_dof[row]) continue;  // 跳过约束DOF
                
                Index row_start = K_csr.row_ptr()[row];
                Index row_end = K_csr.row_ptr()[row + 1];
                
                for (Index k = row_start; k < row_end; ++k) {
                    if (K_csr.col_indices()[k] == i) {
                        // F(row) -= K(row, i) * bc_value
                        F_[row] -= K_csr.values()[k] * bc_val;
                        break;
                    }
                }
            }
        }
    }

    // 2. 修改约束DOF的行和列
    for (Index i = 0; i < n_dofs_; ++i) {
        if (is_bc_dof[i]) {
            // 修改第 i 行: K(i,i)=1, K(i,j≠i)=0
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

            // 修改第 i 列: K(j,i)=0 (对所有 j≠i)
            for (Index row = 0; row < n_dofs_; ++row) {
                if (row == i) continue;
                
                Index start = K_csr.row_ptr()[row];
                Index end = K_csr.row_ptr()[row + 1];
                
                for (Index k = start; k < end; ++k) {
                    if (K_csr.col_indices()[k] == i) {
                        K_csr.values()[k] = 0.0;
                        break;
                    }
                }
            }

            // 设置右端项
            F_[i] = bc_values[i];
        }
    }

    // 转换回 COO 格式
    K_coo_ = csr_to_coo(K_csr);

    FEM_INFO("Applied " + std::to_string(bcs.size()) + " Dirichlet BCs");
}

void Assembler::apply_neumann(const std::vector<NeumannBC>& bcs) {
    if (!assembled_) {
        throw std::runtime_error("Must call assemble() before apply_neumann()");
    }

    std::size_t total_bc_dofs = 0;

    for (const auto& bc : bcs) {
        for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
            const Mesh& mesh = model_.mesh(mesh_id);

            if (!mesh.has_boundary(bc.boundary_name)) {
                continue;
            }

            const auto& boundary_nodes = mesh.boundary(bc.boundary_name);

            if (boundary_nodes.empty()) continue;

            // 当前实现对边界节点进行线积分
            // 注意：这里假设 boundary_nodes 中的节点按空间顺序排列
            // 这依赖于 MeshGenerator::identify_boundaries_2d() 的实现
            // 在矩形域的边界上，节点通常按坐标顺序排列，满足此假设
            
            for (std::size_t i = 1; i < boundary_nodes.size(); ++i) {
                // 处理相邻节点对，模拟边界单元 (Edge2)
                Index n0 = boundary_nodes[i-1];
                Index n1 = boundary_nodes[i];
                
                const auto& c0 = mesh.node(n0).coords();
                const auto& c1 = mesh.node(n1).coords();
                
                // 计算边界段长度
                Real dx = c1[0] - c0[0];
                Real dy = c1[1] - c0[1];
                Real dz = c1[2] - c0[2];
                Real seg_length = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                // 梯形积分：每个节点获得 half-load
                Real half_load = bc.value * seg_length / 2.0;
                
                Index dof0 = n0 * dofs_per_node_ + bc.dof;
                Index dof1 = n1 * dofs_per_node_ + bc.dof;
                
                // 跳过已被 Dirichlet BC 约束的 DOF
                if (dof0 < n_dofs_ && !is_dirichlet_dof_[dof0]) {
                    F_[dof0] += half_load;
                    total_bc_dofs++;
                }
                if (dof1 < n_dofs_ && !is_dirichlet_dof_[dof1]) {
                    F_[dof1] += half_load;
                    total_bc_dofs++;
                }
            }
        }
    }

    FEM_INFO("Applied " + std::to_string(bcs.size()) + " Neumann BCs (" + 
             std::to_string(total_bc_dofs) + " boundary DOFs modified)");
}

SparseMatrixCSR Assembler::matrix() const {
    if (!assembled_) {
        throw std::runtime_error("Must call assemble() before matrix()");
    }
    return coo_to_csr(K_coo_);
}

Vector Assembler::compute_reaction_forces(const std::vector<Real>& u) const {
    if (!has_original_) {
        throw std::runtime_error("Original matrix not saved. apply_dirichlet() not called?");
    }
    
    if (u.size() != n_dofs_) {
        throw std::runtime_error("Displacement vector size mismatch");
    }
    
    // R = K_original * u
    Vector R(n_dofs_, 0.0);
    
    for (Index i = 0; i < n_dofs_; ++i) {
        Index row_start = K_original_.row_ptr()[i];
        Index row_end = K_original_.row_ptr()[i + 1];
        
        Real sum = 0.0;
        for (Index k = row_start; k < row_end; ++k) {
            Index j = K_original_.col_indices()[k];
            sum += K_original_.values()[k] * u[j];
        }
        
        R[i] = sum;
    }
    
    return R;
}

void Assembler::clear() {
    K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);
    F_ = Vector(n_dofs_, 0.0);
    is_dirichlet_dof_.assign(n_dofs_, false);
    assembled_ = false;
    has_original_ = false;
}

} // namespace fem
