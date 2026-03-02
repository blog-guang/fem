/**
 * static_solver.cpp - Static Structural Solver Implementation
 */

#include "assembly/static_solver.h"
#include "core/logger.h"
#include "core/timer.h"
#include <stdexcept>

namespace fem {

StaticSolver::StaticSolver(const Model& model,
                           physics::PhysicsBase* physics,
                           int dofs_per_node)
    : model_(model), physics_(physics), dofs_per_node_(dofs_per_node),
      assembler_(model, dofs_per_node),
      nonlinear_(false)
{
    // 默认：PCG + AMG 预条件
    set_linear_solver("PCG", {{"tol", 1e-8}}, {{"precond", "amg"}});
    
    // 默认非线性求解器参数
    nr_params_.max_iter = 50;
    nr_params_.tol = 1e-6;
    nr_params_.tol_relative = 1e-4;
    nr_params_.line_search_alpha = 1.0;
    nr_params_.verbose = true;
}

void StaticSolver::add_load_step(const LoadStep& step) {
    load_steps_.add_load_step(step);
}

void StaticSolver::clear_load_steps() {
    load_steps_.clear();
}

void StaticSolver::set_linear_solver(
    const std::string& solver_type,
    const SolverFactory::Parameters& params,
    const SolverFactory::StringParams& str_params)
{
    linear_solver_ = SolverFactory::create(solver_type, params, str_params);
}

void StaticSolver::set_nonlinear_params(const NewtonRaphsonParams& params) {
    nr_params_ = params;
}

bool StaticSolver::solve() {
    const Mesh& mesh = model_.mesh(0);  // 使用第一个网格
    
    FEM_INFO("========== Static Structural Solver ==========");
    FEM_INFO("Model: " + std::to_string(mesh.num_nodes()) + " nodes, " +
             std::to_string(mesh.num_elements()) + " elements");
    FEM_INFO("DOFs per node: " + std::to_string(dofs_per_node_));
    FEM_INFO("Total DOFs: " + std::to_string(mesh.num_nodes() * dofs_per_node_));
    FEM_INFO("Load steps: " + std::to_string(load_steps_.num_load_steps()));
    FEM_INFO("Nonlinear: " + std::string(nonlinear_ ? "YES" : "NO"));
    FEM_INFO("==============================================");
    
    // 初始化解向量（如果需要）
    std::size_t total_dofs = mesh.num_nodes() * dofs_per_node_;
    if (u_.size() != total_dofs) {
        u_ = Vector(total_dofs, 0.0);
    }
    
    // 清空历史
    history_.clear();
    
    // 循环所有载荷步
    bool all_converged = true;
    for (std::size_t i = 0; i < load_steps_.num_load_steps(); i++) {
        const LoadStep& step = load_steps_.load_steps()[i];
        
        FEM_INFO("--------------------------------------------------");
        FEM_INFO("Solving LoadStep " + std::to_string(step.id()));
        step.print();
        
        bool converged = solve_load_step(step.id());
        
        if (!converged) {
            all_converged = false;
            FEM_ERROR("LoadStep " + std::to_string(step.id()) + " FAILED");
            break;  // 停止求解
        }
    }
    
    // 打印摘要
    FEM_INFO("==============================================");
    print_summary();
    
    return all_converged;
}

bool StaticSolver::solve_load_step(int load_step_id) {
    const LoadStep& step = load_steps_.get_load_step(load_step_id);
    
    // 循环所有子步
    bool all_converged = true;
    for (int substep = 1; substep <= step.num_substeps(); substep++) {
        bool converged = solve_substep(step, substep);
        
        if (!converged) {
            all_converged = false;
            FEM_ERROR("Substep " + std::to_string(substep) + " FAILED");
            break;
        }
    }
    
    return all_converged;
}

bool StaticSolver::solve_substep(const LoadStep& step, int substep_id) {
    Timer timer;
    timer.start();
    
    // 计算当前时间
    Real time = step.substep_time(substep_id);
    
    // 获取边界条件和载荷（插值到当前时间）
    auto bcs = step.get_displacements(time);
    auto forces = step.get_forces(time);
    
    // 选择求解方法
    bool converged = false;
    if (nonlinear_) {
        converged = solve_nonlinear(bcs, forces);
    } else {
        converged = solve_linear(bcs, forces);
    }
    
    timer.stop();
    
    // 记录收敛历史
    SubstepResult result;
    result.load_step_id = step.id();
    result.substep_id = substep_id;
    result.time = time;
    result.converged = converged;
    result.iterations = 1;  // 线性求解：1 次迭代
    result.residual_norm = 0.0;  // TODO: 提取残差
    result.displacement_norm = u_.norm();
    result.solve_time = timer.elapsed_s();
    
    history_.add_substep(result);
    
    return converged;
}

bool StaticSolver::solve_linear(
    const std::map<std::pair<Index, int>, Real>& bcs,
    const std::map<std::pair<Index, int>, Real>& forces)
{
    // 装配刚度矩阵和载荷向量
    assembler_.assemble([this](Index elem_id, const Mesh& mesh,
                               DenseMatrix& Ke, Vector& Fe) {
        physics_->compute_element(elem_id, mesh, Ke, Fe);
    });
    
    // 施加力载荷
    for (const auto& [key, value] : forces) {
        Index node_id = key.first;
        int dof = key.second;
        Index global_dof = node_id * dofs_per_node_ + dof;
        assembler_.rhs()[global_dof] += value;
    }
    
    // 施加位移边界条件
    for (const auto& [key, value] : bcs) {
        Index node_id = key.first;
        int dof = key.second;
        Index global_dof = node_id * dofs_per_node_ + dof;
        assembler_.apply_dirichlet_single(global_dof, value);
    }
    
    // 求解 Ku = F
    SparseMatrixCSR K = assembler_.matrix();
    Vector& F = assembler_.rhs();
    
    auto result = linear_solver_->solve(K, F, u_);
    
    return result.converged;
}

bool StaticSolver::solve_nonlinear(
    const std::map<std::pair<Index, int>, Real>& bcs,
    const std::map<std::pair<Index, int>, Real>& forces)
{
    // TODO: 实现完整的 Newton-Raphson 非线性求解
    // 当前暂时使用线性求解
    FEM_WARN("Nonlinear solve not fully implemented yet, using linear solve");
    return solve_linear(bcs, forces);
}

void StaticSolver::print_summary() const {
    history_.print_summary();
}

void StaticSolver::export_history(const std::string& filename) const {
    history_.export_csv(filename);
}

}  // namespace fem
