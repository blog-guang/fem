/**
 * pcg.h - Preconditioned Conjugate Gradient Solver
 * 
 * 预条件共轭梯度法求解器
 * 
 * 算法：PCG with various preconditioners
 * 
 * 预条件器选项：
 * - Jacobi (对角预条件器): M = diag(K)
 * - SSOR (对称超松弛): M = (D+L)D^{-1}(D+U)
 * - ILU(0) (不完全 LU 分解)
 * - AMG (代数多重网格，使用 AMGCL 库)
 * 
 * 默认使用 Jacobi 预条件器（最简单且稳定）
 */

#pragma once

#include "math/solver.h"
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
 * ILU(0) (不完全 LU 分解) 预条件器
 * 
 * A ≈ L * U (保持原始稀疏模式，不产生新的非零元)
 * 
 * 算法：
 * 1. 复制 A 的非零模式到 LU
 * 2. 对 LU 进行原位 ILU(0) 分解
 * 3. 应用时求解 L*U*z = r (前向+后向替换)
 * 
 * 优点：
 * - 收敛速度快（通常比 Jacobi 快 2-5 倍）
 * - 适合大型稀疏系统
 * - 内存开销小（仅存储原矩阵非零元）
 * 
 * 缺点：
 * - 构建时间较长
 * - 串行算法
 * - 可能不稳定（需要对角占优）
 */
class ILUPreconditioner : public Preconditioner {
public:
    void apply(const std::vector<Real>& r, std::vector<Real>& z) const override;
    void build(const SparseMatrixCSR& K) override;

private:
    SparseMatrixCSR LU_;  // 存储 L 和 U (L 的对角线隐式为 1)
    
    /**
     * 前向替换: 求解 L*y = r
     * (L 的对角线为 1)
     */
    void forward_solve(const std::vector<Real>& r, std::vector<Real>& y) const;
    
    /**
     * 后向替换: 求解 U*z = y
     */
    void backward_solve(const std::vector<Real>& y, std::vector<Real>& z) const;
};

/**
 * AMG (代数多重网格) 预条件器
 * 
 * 使用 AMGCL 库实现
 * 
 * 优点：
 * - 收敛速度非常快（O(N) 复杂度）
 * - 适合大规模问题（百万级自由度）
 * - 对网格质量不敏感
 * - 自动构建多层网格
 * 
 * 缺点：
 * - 构建时间较长
 * - 内存开销较大
 * - 需要 AMGCL 库
 */
class AMGPreconditioner : public Preconditioner {
public:
    AMGPreconditioner();
    ~AMGPreconditioner();
    
    void apply(const std::vector<Real>& r, std::vector<Real>& z) const override;
    void build(const SparseMatrixCSR& K) override;

private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;  // PIMPL idiom to hide AMGCL headers
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
     * @param precond_type 预条件器类型 ("jacobi", "ssor", "ilu", "amg", "none")
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
