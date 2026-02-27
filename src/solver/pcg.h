/**
 * pcg.h - Preconditioned Conjugate Gradient Solver
 * 
 * 预条件共轭梯度法求解器
 * 
 * 算法：PCG with Jacobi (diagonal) preconditioner
 * 
 * 预条件器选项：
 * - Jacobi (对角预条件器): M = diag(K)
 * - SSOR (对称超松弛): M = (D+L)D^{-1}(D+U)
 * - ILU(0) (不完全 LU 分解)
 * 
 * 默认使用 Jacobi 预条件器（最简单且稳定）
 */

#pragma once

#include "solver/solver.h"
#include <memory>

namespace fem {

/**
 * 预条件器基类
 */
class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    
    /**
     * 应用预条件器: z = M^{-1} * r
     */
    virtual void apply(const std::vector<Real>& r, std::vector<Real>& z) const = 0;
    
    /**
     * 构建预条件器
     */
    virtual void build(const SparseMatrixCSR& K) = 0;
};

/**
 * Jacobi (对角) 预条件器
 * 
 * M = diag(K)
 * M^{-1} * r = r ./ diag(K)
 * 
 * 优点：简单、并行、内存少
 * 缺点：收敛速度一般
 */
class JacobiPreconditioner : public Preconditioner {
public:
    void apply(const std::vector<Real>& r, std::vector<Real>& z) const override;
    void build(const SparseMatrixCSR& K) override;

private:
    std::vector<Real> diag_inv_;  // 对角元素的倒数
};

/**
 * SSOR (对称超松弛) 预条件器
 * 
 * M = (D + ωL) D^{-1} (D + ωU)
 * 
 * 其中：
 * - D: 对角矩阵
 * - L: 严格下三角
 * - U: 严格上三角
 * - ω: 松弛因子 (通常 ω = 1)
 * 
 * 优点：收敛速度较快
 * 缺点：串行算法、实现复杂
 */
class SSORPreconditioner : public Preconditioner {
public:
    explicit SSORPreconditioner(Real omega = 1.0) : omega_(omega) {}
    
    void apply(const std::vector<Real>& r, std::vector<Real>& z) const override;
    void build(const SparseMatrixCSR& K) override;

private:
    Real omega_;                    // 松弛因子
    const SparseMatrixCSR* K_;      // 指向原始矩阵
    std::vector<Real> diag_inv_;    // 对角元素的倒数
};

/**
 * 预条件共轭梯度求解器 (PCG)
 * 
 * 算法：
 * 
 * 1. r0 = b - A*x0
 * 2. z0 = M^{-1} * r0
 * 3. p0 = z0
 * 4. for k = 0, 1, 2, ...
 *      α_k = (r_k, z_k) / (p_k, A*p_k)
 *      x_{k+1} = x_k + α_k * p_k
 *      r_{k+1} = r_k - α_k * A*p_k
 *      z_{k+1} = M^{-1} * r_{k+1}
 *      β_k = (r_{k+1}, z_{k+1}) / (r_k, z_k)
 *      p_{k+1} = z_{k+1} + β_k * p_k
 * 
 * 相比标准 CG:
 * - 使用 z = M^{-1} * r 代替 r
 * - 内积从 (r,r) 变为 (r,z)
 */
class PCGSolver : public LinearSolver {
public:
    /**
     * 构造函数
     * 
     * @param precond_type 预条件器类型 ("jacobi", "ssor", "none")
     * @param omega SSOR 松弛因子 (仅当 precond_type = "ssor" 时有效)
     */
    explicit PCGSolver(const std::string& precond_type = "jacobi", Real omega = 1.0);
    
    SolveResult solve(const SparseMatrixCSR& K,
                     const std::vector<Real>& F,
                     std::vector<Real>& x) override;

private:
    std::string precond_type_;
    Real omega_;
    std::unique_ptr<Preconditioner> precond_;
};

}  // namespace fem
