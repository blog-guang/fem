/**
 * nonlinear_solver.cpp - 非线性求解器实现
 */

#include "solver/nonlinear_solver.h"
#include "math/solver_factory.h"
#include "core/logger.h"
#include "core/timer.h"
#include <cmath>
#include <sstream>
#include <iomanip>

namespace fem {
namespace solver {

NonlinearSolver::NonlinearSolver(const Model& model, Index dofs_per_node)
    : model_(model)
    , dofs_per_node_(dofs_per_node)
    , physics_(nullptr)
    , verbose_(false)
    , load_scale_(1.0)
{
    // 计算总自由度数
    total_dofs_ = 0;
    for (int i = 0; i < model_.num_meshes(); i++) {
        total_dofs_ += model_.mesh(i).num_nodes() * dofs_per_node_;
    }
    
    // 默认使用 CG 求解器
    linear_solver_ = SolverFactory::create("CG", {{"tolerance", 1e-8}});
}

// ═══════════════════════════════════════════════════════════
// 边界条件
// ═══════════════════════════════════════════════════════════

void NonlinearSolver::apply_dirichlet(Index node_id, Index dof, Real value) {
    Index global_dof = node_id * dofs_per_node_ + dof;
    
    if (global_dof < total_dofs_) {
        dirichlet_bc_[global_dof] = value;
    }
}

void NonlinearSolver::apply_dirichlet(const std::vector<DirichletBC>& bcs) {
    for (const auto& bc : bcs) {
        // TODO: 处理命名边界
        // 当前仅支持直接指定节点 ID 和 DOF
    }
}

// ═══════════════════════════════════════════════════════════
// 求解
// ═══════════════════════════════════════════════════════════

bool NonlinearSolver::solve(const Vector& F_ext, const Vector& u_initial) {
    if (!physics_) {
        FEM_ERROR("Physics module not set!");
        return false;
    }
    
    // 初始化
    F_ext_ = F_ext * load_scale_;
    
    if (u_initial.size() == F_ext.size()) {
        u_current_ = u_initial;
    } else {
        u_current_.resize(F_ext.size(), 0.0);
    }
    
    // 应用 Dirichlet 边界条件到初始位移
    for (const auto& [dof, value] : dirichlet_bc_) {
        if (dof < u_current_.size()) {
            u_current_[dof] = value;
        }
    }
    
    history_.clear();
    
    if (verbose_) {
        std::ostringstream oss;
        oss << "=== Nonlinear Solver (Newton-Raphson) ===\n"
            << "  Total DOFs: " << F_ext.size() << "\n"
            << "  Prescribed DOFs: " << dirichlet_bc_.size() << "\n"
            << "  Load scale: " << load_scale_ << "\n"
            << "  Max iterations: " << criteria_.max_iterations << "\n"
            << "  Residual tol: " << std::scientific << criteria_.residual_tol << "\n"
            << "  Displacement tol: " << criteria_.displacement_tol;
        FEM_INFO(oss.str());
    }
    
    // Newton-Raphson 迭代
    for (int iter = 0; iter < criteria_.max_iterations; iter++) {
        bool converged = iterate();
        
        if (converged) {
            if (verbose_) {
                std::ostringstream oss;
                oss << "Newton-Raphson converged in " << iter + 1 << " iterations";
                FEM_INFO(oss.str());
            }
            return true;
        }
    }
    
    if (verbose_) {
        std::ostringstream oss;
        oss << "Newton-Raphson did not converge in " << criteria_.max_iterations << " iterations";
        FEM_WARN(oss.str());
    }
    
    return false;
}

bool NonlinearSolver::iterate() {
    NonlinearIteration iter;
    iter.iteration = static_cast<int>(history_.size());
    
    Timer timer;
    
    // 1. 计算残差：R = F_int - F_ext
    Vector R;
    compute_residual(u_current_, R);
    
    // 2. 计算切线刚度：K_t
    SparseMatrixCSR K_t;
    compute_tangent_stiffness(u_current_, K_t);
    
    // 3. 应用边界条件
    apply_bc_to_system(K_t, R);
    
    // 4. 求解：K_t * du = -R
    Vector du;
    Vector neg_R = R * (-1.0);
    
    linear_solver_->solve(K_t, neg_R, du);
    
    iter.solve_time = timer.elapsed_s();
    
    // 5. 更新位移：u^(k+1) = u^k + du
    u_current_ = u_current_ + du;
    
    // 重新应用 Dirichlet 边界条件（确保固定）
    for (const auto& [dof, value] : dirichlet_bc_) {
        if (dof < u_current_.size()) {
            u_current_[dof] = value;
        }
    }
    
    // 6. 检查收敛
    bool converged = check_convergence(R, du, iter);
    
    history_.push_back(iter);
    
    if (verbose_) {
        print_iteration(iter);
    }
    
    return converged;
}

void NonlinearSolver::reset() {
    history_.clear();
    u_current_.clear();
    F_ext_.clear();
    dirichlet_bc_.clear();
    load_scale_ = 1.0;
}

// ═══════════════════════════════════════════════════════════
// 内部计算
// ═══════════════════════════════════════════════════════════

void NonlinearSolver::compute_residual(const Vector& u, Vector& R) {
    // R = F_int(u) - F_ext
    Vector F_int;
    compute_internal_force(u, F_int);
    
    R = F_int - F_ext_;
}

void NonlinearSolver::compute_internal_force(const Vector& u, Vector& F_int) {
    // 初始化
    F_int.resize(u.size(), 0.0);
    
    // 遍历所有网格
    for (int mesh_id = 0; mesh_id < model_.num_meshes(); mesh_id++) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        // 遍历所有单元
        for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
            // 计算单元内力
            Vector F_elem = physics_->compute_internal_force(elem_id, mesh, u);
            
            // 装配到全局向量
            const auto& elem = mesh.element(elem_id);
            const auto& node_ids = elem.nodes();
            
            for (std::size_t i = 0; i < node_ids.size(); i++) {
                Index node_id = node_ids[i];
                
                for (Index d = 0; d < dofs_per_node_; d++) {
                    Index global_dof = node_id * dofs_per_node_ + d;
                    Index local_dof = i * dofs_per_node_ + d;
                    
                    if (global_dof < F_int.size() && local_dof < F_elem.size()) {
                        F_int[global_dof] += F_elem[local_dof];
                    }
                }
            }
        }
    }
}

void NonlinearSolver::compute_tangent_stiffness(
    const Vector& u,
    SparseMatrixCSR& K_t)
{
    // 使用 Assembler + 回调函数方式
    Assembler assembler(model_, dofs_per_node_);
    
    // 遍历所有网格并装配
    for (int mesh_id = 0; mesh_id < model_.num_meshes(); mesh_id++) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        // 创建回调函数（捕获当前位移和物理模块）
        auto element_callback = [this, &mesh, &u](
            Index elem_id,
            const Mesh& /* unused_mesh */,
            DenseMatrix& Ke,
            Vector& Fe) {
            
            // 调用几何非线性物理模块
            physics_->compute_element_nonlinear(elem_id, mesh, u, Ke, Fe);
        };
        
        // 使用回调装配刚度矩阵
        assembler.assemble_matrix(element_callback);
    }
    
    // 获取刚度矩阵（CSR 格式）
    K_t = assembler.matrix();
}

void NonlinearSolver::apply_bc_to_system(SparseMatrixCSR& K, Vector& R) {
    // 应用 Dirichlet 边界条件
    // 方法：置1法
    // - K(i,i) = 1
    // - K(i,j) = 0 (j ≠ i)
    // - R(i) = 0 (位移已知)
    
    for (const auto& [dof, value] : dirichlet_bc_) {
        if (dof < R.size()) {
            R[dof] = 0.0;  // 残差置零
        }
        
        // 修改刚度矩阵（简化：仅置零残差，依赖求解器稳定性）
        // 完整实现需要修改 CSR 矩阵结构
    }
}

bool NonlinearSolver::check_convergence(
    const Vector& R,
    const Vector& du,
    NonlinearIteration& iter)
{
    // 计算残差范数
    iter.residual_norm = 0.0;
    for (std::size_t i = 0; i < R.size(); i++) {
        iter.residual_norm += R[i] * R[i];
    }
    iter.residual_norm = std::sqrt(iter.residual_norm);
    
    // 计算位移增量范数
    iter.displacement_norm = 0.0;
    for (std::size_t i = 0; i < du.size(); i++) {
        iter.displacement_norm += du[i] * du[i];
    }
    iter.displacement_norm = std::sqrt(iter.displacement_norm);
    
    // 计算能量：|ΔuᵀR|
    iter.energy = 0.0;
    for (std::size_t i = 0; i < std::min(du.size(), R.size()); i++) {
        iter.energy += du[i] * R[i];
    }
    iter.energy = std::abs(iter.energy);
    
    // 检查收敛准则
    bool converged = true;
    
    if (criteria_.use_residual) {
        if (iter.residual_norm > criteria_.residual_tol) {
            converged = false;
        }
    }
    
    if (criteria_.use_displacement) {
        if (iter.displacement_norm > criteria_.displacement_tol) {
            converged = false;
        }
    }
    
    if (criteria_.use_energy) {
        if (iter.energy > criteria_.energy_tol) {
            converged = false;
        }
    }
    
    iter.converged = converged;
    return converged;
}

void NonlinearSolver::print_iteration(const NonlinearIteration& iter) {
    std::ostringstream oss;
    oss << "  Iter " << std::setw(3) << iter.iteration
        << ": |R| = " << std::scientific << std::setprecision(3) << iter.residual_norm
        << ", |du| = " << iter.displacement_norm
        << ", E = " << iter.energy
        << " (" << std::fixed << std::setprecision(3) << iter.solve_time << "s)";
    
    if (iter.converged) {
        oss << " ✓";
    }
    
    FEM_INFO(oss.str());
}

}  // namespace solver
}  // namespace fem
