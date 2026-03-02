/**
 * newton_raphson_solver.cpp - Newton-Raphson 非线性求解器实现
 */

#include "solver/newton_raphson_solver.h"
#include "math/solver_factory.h"
#include "physics/elasticity_nonlinear.h"
#include "core/logger.h"
#include <cmath>
#include <iomanip>
#include <sstream>

namespace fem {
namespace solver {

NewtonRaphsonSolver::NewtonRaphsonSolver(
    Model& model,
    physics::PhysicsBase* physics,
    int dofs_per_node)
    : model_(model)
    , physics_(physics)
    , dofs_per_node_(dofs_per_node)
    , verbose_(false)
    , load_scale_(1.0)
{
    // 默认使用 CG 求解器
    linear_solver_ = SolverFactory::create("CG", {{"tolerance", 1e-8}});
    
    // 计算总节点数
    Index total_nodes = 0;
    for (int i = 0; i < model_.num_meshes(); i++) {
        total_nodes += model_.mesh(i).num_nodes();
    }
    
    // 创建装配器
    assembler_ = std::make_unique<Assembler>(
        total_nodes,
        dofs_per_node_
    );
}

// ═══════════════════════════════════════════════════════════
// 求解控制
// ═══════════════════════════════════════════════════════════

bool NewtonRaphsonSolver::solve(const Vector& F_ext, const Vector& u_initial) {
    // 初始化
    F_ext_ = F_ext * load_scale_;
    
    if (u_initial.size() == F_ext.size()) {
        u_current_ = u_initial;
    } else {
        u_current_.resize(F_ext.size(), 0.0);
    }
    
    history_.clear();
    
    if (verbose_) {
        std::ostringstream oss;
        oss << "Newton-Raphson solver started\n"
            << "  DOFs: " << F_ext.size() << "\n"
            << "  Load scale: " << load_scale_ << "\n"
            << "  Max iterations: " << criteria_.max_iterations << "\n"
            << "  Residual tol: " << criteria_.residual_tol << "\n"
            << "  Displacement tol: " << criteria_.displacement_tol;
        LOG_INFO(oss.str());
    }
    
    // Newton-Raphson 迭代
    for (int iter = 0; iter < criteria_.max_iterations; iter++) {
        bool converged = iterate();
        
        if (converged) {
            if (verbose_) {
                std::ostringstream oss;
                oss << "Newton-Raphson converged in " << iter + 1 << " iterations";
                LOG_INFO(oss.str());
            }
            return true;
        }
    }
    
    if (verbose_) {
        std::ostringstream oss;
        oss << "Newton-Raphson did not converge in " << criteria_.max_iterations << " iterations";
        LOG_WARNING(oss.str());
    }
    
    return false;
}

bool NewtonRaphsonSolver::iterate() {
    NewtonRaphsonIteration iter;
    iter.iteration = history_.size();
    
    // 1. 计算残差：R = F_int - F_ext
    Vector R;
    compute_residual(u_current_, R);
    
    // 2. 计算切线刚度：K_t
    SparseMatrix K_t;
    compute_tangent_stiffness(u_current_, K_t);
    
    // 3. 应用边界条件
    apply_boundary_conditions(R, K_t);
    
    // 4. 求解：K_t * du = -R
    Vector du;
    Vector neg_R = R * (-1.0);
    
    linear_solver_->solve(K_t, neg_R, du);
    
    // 5. 更新位移：u^(k+1) = u^k + du
    u_current_ = u_current_ + du;
    
    // 6. 检查收敛
    bool converged = check_convergence(R, du, iter);
    
    history_.push_back(iter);
    
    if (verbose_) {
        print_iteration(iter);
    }
    
    return converged;
}

// ═══════════════════════════════════════════════════════════
// 内部计算
// ═══════════════════════════════════════════════════════════

void NewtonRaphsonSolver::compute_residual(const Vector& u, Vector& R) {
    // R = F_int(u) - F_ext
    Vector F_int;
    compute_internal_force(u, F_int);
    
    R = F_int - F_ext_;
}

void NewtonRaphsonSolver::compute_internal_force(const Vector& u, Vector& F_int) {
    // 初始化
    F_int.resize(u.size(), 0.0);
    
    // 检查物理模块是否支持几何非线性
    auto* nonlinear_physics = dynamic_cast<physics::ElasticityNonlinear*>(physics_);
    
    if (!nonlinear_physics) {
        throw std::runtime_error(
            "NewtonRaphsonSolver requires ElasticityNonlinear physics module");
    }
    
    // 遍历所有网格
    for (int mesh_id = 0; mesh_id < model_.num_meshes(); mesh_id++) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        // 遍历所有单元
        for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
            // 计算单元内力
            Vector F_elem = nonlinear_physics->compute_internal_force(
                elem_id, mesh, u);
            
            // 装配到全局向量
            const auto& elem = mesh.element(elem_id);
            const auto& node_ids = elem.nodes();
            
            for (std::size_t i = 0; i < node_ids.size(); i++) {
                Index node_id = node_ids[i];
                
                for (int d = 0; d < dofs_per_node_; d++) {
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

void NewtonRaphsonSolver::compute_tangent_stiffness(
    const Vector& u,
    SparseMatrixCSR& K_t)
{
    // 检查物理模块
    auto* nonlinear_physics = dynamic_cast<physics::ElasticityNonlinear*>(physics_);
    
    if (!nonlinear_physics) {
        throw std::runtime_error(
            "NewtonRaphsonSolver requires ElasticityNonlinear physics module");
    }
    
    // 重置装配器
    assembler_->reset();
    
    // 遍历所有网格
    for (int mesh_id = 0; mesh_id < model_.num_meshes(); mesh_id++) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        // 遍历所有单元
        for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
            // 计算单元切线刚度（包含几何刚度）
            DenseMatrix Ke;
            Vector Fe;  // 忽略（已在 F_int 中计算）
            
            nonlinear_physics->compute_element_nonlinear(
                elem_id, mesh, u, Ke, Fe);
            
            // 装配到全局刚度矩阵
            assembler_->add_element_matrix(Ke, mesh.element(elem_id));
        }
    }
    
    // 获取全局刚度矩阵（转换为 CSR 格式）
    K_t = assembler_->stiffness_csr();
}

void NewtonRaphsonSolver::apply_boundary_conditions(
    Vector& R,
    SparseMatrixCSR& K_t)
{
    // 应用边界条件到刚度矩阵和残差向量
    // 对于位移边界条件：
    // - 刚度矩阵对角元设为1
    // - 残差向量对应位置设为0（位移已知）
    
    const auto& bc_map = assembler_->dirichlet_bc();
    
    for (const auto& [global_dof, value] : bc_map) {
        if (global_dof < R.size()) {
            // 残差设为0（位移已固定）
            R[global_dof] = 0.0;
            
            // 刚度矩阵对应行列处理
            // 简化处理：对角元设为1，其他元设为0
            // 完整实现需要修改稀疏矩阵结构
            // TODO: 实现 CSR 矩阵的边界条件施加
        }
    }
}

bool NewtonRaphsonSolver::check_convergence(
    const Vector& R,
    const Vector& du,
    NewtonRaphsonIteration& iter)
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
    
    // 计算能量：ΔuᵀR
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

void NewtonRaphsonSolver::print_iteration(const NewtonRaphsonIteration& iter) {
    std::ostringstream oss;
    oss << "  Iter " << std::setw(3) << iter.iteration
        << ": |R| = " << std::scientific << std::setprecision(3) << iter.residual_norm
        << ", |du| = " << iter.displacement_norm
        << ", E = " << iter.energy;
    
    if (iter.converged) {
        oss << " [CONVERGED]";
    }
    
    LOG_INFO(oss.str());
}

// ═══════════════════════════════════════════════════════════
// 辅助功能
// ═══════════════════════════════════════════════════════════

void NewtonRaphsonSolver::reset() {
    history_.clear();
    u_current_.clear();
    F_ext_.clear();
    load_scale_ = 1.0;
}

}  // namespace solver
}  // namespace fem
