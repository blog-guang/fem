/**
 * newton_raphson.h - Newton-Raphson 非线性求解器
 * 
 * 求解非线性方程组: R(u) = 0
 * 其中 R(u) = F_int(u) - F_ext
 * 
 * 迭代公式:
 *   K_t(u_n) * Δu = -R(u_n)
 *   u_{n+1} = u_n + Δu
 * 
 * 其中 K_t = ∂R/∂u 是切线刚度矩阵 (Tangent stiffness)
 */

#pragma once

#include "core/types.h"
#include "math/sparse_matrix.h"
#include "math/vector.h"
#include "solver/cg.h"

#include <functional>
#include <vector>

namespace fem {

/**
 * 非线性问题接口
 * 
 * 用户需要实现：
 * 1. compute_residual: 计算残差向量 R(u)
 * 2. compute_tangent: 计算切线刚度矩阵 K_t(u)
 */
struct NonlinearProblem {
    /**
     * 计算残差向量 R(u) = F_int(u) - F_ext
     * @param u 当前解向量
     * @param R 输出：残差向量
     */
    virtual void compute_residual(const std::vector<Real>& u, 
                                  std::vector<Real>& R) = 0;
    
    /**
     * 计算切线刚度矩阵 K_t(u) = ∂R/∂u
     * @param u 当前解向量
     * @param K_t 输出：切线刚度矩阵
     */
    virtual void compute_tangent(const std::vector<Real>& u,
                                SparseMatrixCSR& K_t) = 0;
    
    virtual ~NonlinearProblem() = default;
};

/**
 * Newton-Raphson 求解器参数
 */
struct NewtonRaphsonParams {
    Real tol = 1e-6;              ///< 残差收敛容差
    Real tol_relative = 1e-6;     ///< 相对收敛容差
    Index max_iter = 50;          ///< 最大迭代次数
    Real line_search_alpha = 1.0; ///< 线搜索步长（1.0 = 不使用线搜索）
    bool verbose = false;         ///< 是否输出详细信息
};

/**
 * Newton-Raphson 求解结果
 */
struct NewtonRaphsonResult {
    bool converged = false;       ///< 是否收敛
    Index iterations = 0;         ///< 迭代次数
    Real residual_norm = 0.0;     ///< 最终残差范数
    Real initial_residual = 0.0;  ///< 初始残差范数
    std::vector<Real> residual_history;  ///< 残差历史
};

/**
 * Newton-Raphson 非线性求解器
 * 
 * 使用方法：
 * ```cpp
 * NewtonRaphsonSolver solver;
 * solver.set_params(params);
 * 
 * std::vector<Real> u0(n, 0.0);  // 初始猜测
 * auto result = solver.solve(problem, u0);
 * 
 * if (result.converged) {
 *     // u0 now contains the solution
 * }
 * ```
 */
class NewtonRaphsonSolver {
public:
    /**
     * 设置求解参数
     */
    void set_params(const NewtonRaphsonParams& params) {
        params_ = params;
    }
    
    /**
     * 求解非线性问题
     * @param problem 非线性问题接口
     * @param u 输入/输出：初始猜测和最终解
     * @return 求解结果
     */
    NewtonRaphsonResult solve(NonlinearProblem& problem, 
                             std::vector<Real>& u);
    
private:
    NewtonRaphsonParams params_;
    
    /**
     * 计算向量的 L2 范数
     */
    Real compute_norm(const std::vector<Real>& v) const;
    
    /**
     * 线搜索（可选）
     * 寻找步长 α 使得 ||R(u + α*Δu)|| < ||R(u)||
     */
    Real line_search(NonlinearProblem& problem,
                    const std::vector<Real>& u,
                    const std::vector<Real>& du,
                    Real residual_norm);
};

} // namespace fem
