#include "math/newton_raphson.h"
#include "core/logger.h"
#include <cmath>
#include <algorithm>

namespace fem {

NewtonRaphsonResult NewtonRaphsonSolver::solve(NonlinearProblem& problem,
                                                std::vector<Real>& u) {
    NewtonRaphsonResult result;
    
    const Index n = u.size();
    std::vector<Real> R(n);
    std::vector<Real> du(n);
    
    // 计算初始残差
    problem.compute_residual(u, R);
    result.initial_residual = compute_norm(R);
    result.residual_norm = result.initial_residual;
    
    if (params_.verbose) {
        FEM_INFO("Newton-Raphson: Initial residual = " + fmt_sci(result.initial_residual));
    }
    
    // 如果初始残差已经很小，直接返回
    if (result.initial_residual < params_.tol) {
        result.converged = true;
        result.iterations = 0;
        return result;
    }
    
    // Newton-Raphson 迭代
    for (Index iter = 0; iter < params_.max_iter; ++iter) {
        // 1. 计算切线刚度矩阵 K_t(u_n)
        SparseMatrixCSR K_t;
        problem.compute_tangent(u, K_t);
        
        // 2. 求解线性系统: K_t * du = -R
        std::vector<Real> neg_R(n);
        for (Index i = 0; i < n; ++i) {
            neg_R[i] = -R[i];
        }
        
        // 使用 CG 求解器求解线性系统
        CGSolver linear_solver;
        linear_solver.set_tol(1e-8);
        linear_solver.set_max_iter(10000);
        
        auto linear_result = linear_solver.solve(K_t, neg_R, du);
        
        if (!linear_result.converged) {
            FEM_WARN("Newton-Raphson: Linear solver failed at iteration " + 
                    std::to_string(iter));
            result.converged = false;
            result.iterations = iter + 1;
            return result;
        }
        
        // 3. 线搜索（可选）
        Real alpha = params_.line_search_alpha;
        if (alpha < 1.0) {
            alpha = line_search(problem, u, du, result.residual_norm);
        }
        
        // 4. 更新解: u_{n+1} = u_n + α * du
        for (Index i = 0; i < n; ++i) {
            u[i] += alpha * du[i];
        }
        
        // 5. 计算新的残差
        problem.compute_residual(u, R);
        result.residual_norm = compute_norm(R);
        result.residual_history.push_back(result.residual_norm);
        
        if (params_.verbose) {
            FEM_INFO("  Iter " + std::to_string(iter + 1) + 
                    ": residual = " + fmt_sci(result.residual_norm) +
                    ", alpha = " + std::to_string(alpha));
        }
        
        // 6. 检查收敛
        bool converged_abs = (result.residual_norm < params_.tol);
        bool converged_rel = (result.residual_norm < params_.tol_relative * result.initial_residual);
        
        if (converged_abs || converged_rel) {
            result.converged = true;
            result.iterations = iter + 1;
            
            if (params_.verbose) {
                FEM_INFO("Newton-Raphson converged in " + std::to_string(result.iterations) + 
                        " iterations");
            }
            
            return result;
        }
    }
    
    // 未收敛
    result.converged = false;
    result.iterations = params_.max_iter;
    
    FEM_WARN("Newton-Raphson failed to converge after " + 
            std::to_string(params_.max_iter) + " iterations");
    FEM_WARN("  Final residual: " + fmt_sci(result.residual_norm));
    
    return result;
}

Real NewtonRaphsonSolver::compute_norm(const std::vector<Real>& v) const {
    Real sum = 0.0;
    for (Real x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

Real NewtonRaphsonSolver::line_search(NonlinearProblem& problem,
                                      const std::vector<Real>& u,
                                      const std::vector<Real>& du,
                                      Real residual_norm) {
    // 简单的回溯线搜索（Backtracking line search）
    constexpr Real beta = 0.5;    // 收缩因子
    constexpr Real c = 1e-4;      // Armijo 条件参数
    constexpr Index max_iter = 10;
    
    Real alpha = 1.0;
    std::vector<Real> u_new(u.size());
    std::vector<Real> R_new(u.size());
    
    for (Index iter = 0; iter < max_iter; ++iter) {
        // u_new = u + α * du
        for (std::size_t i = 0; i < u.size(); ++i) {
            u_new[i] = u[i] + alpha * du[i];
        }
        
        // 计算新残差
        problem.compute_residual(u_new, R_new);
        Real residual_new = compute_norm(R_new);
        
        // Armijo 条件：||R_new|| < (1 - c*α) * ||R||
        if (residual_new < (1.0 - c * alpha) * residual_norm) {
            return alpha;
        }
        
        // 收缩步长
        alpha *= beta;
    }
    
    // 如果线搜索失败，返回最小步长
    return alpha;
}

} // namespace fem
