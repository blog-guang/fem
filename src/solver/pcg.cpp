/**
 * pcg.cpp - PCG Solver Implementation
 */

#include "solver/pcg.h"
#include "core/logger.h"
#include <cmath>
#include <algorithm>

namespace fem {

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

static Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

static void axpy(Real alpha, const std::vector<Real>& x, std::vector<Real>& y) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] += alpha * x[i];
    }
}

// ═══════════════════════════════════════════════════════════
// Jacobi 预条件器
// ═══════════════════════════════════════════════════════════

void JacobiPreconditioner::build(const SparseMatrixCSR& K) {
    std::size_t n = K.rows();
    diag_inv_.resize(n);
    
    // 提取对角元素并取倒数
    for (std::size_t i = 0; i < n; ++i) {
        Real diag = 0.0;
        
        // 找到对角元素 K(i,i)
        for (std::size_t idx = K.row_ptr()[i]; idx < K.row_ptr()[i + 1]; ++idx) {
            if (K.col_indices()[idx] == i) {
                diag = K.values()[idx];
                break;
            }
        }
        
        // 检查对角元素是否为零或太小
        if (std::abs(diag) < 1e-30) {
            FEM_WARN("Jacobi preconditioner: diagonal element " + std::to_string(i) +
                    " is near zero (" + std::to_string(diag) + "), using 1.0");
            diag_inv_[i] = 1.0;
        } else {
            diag_inv_[i] = 1.0 / diag;
        }
    }
    
    FEM_INFO("Jacobi preconditioner built: " + std::to_string(n) + " DOFs");
}

void JacobiPreconditioner::apply(const std::vector<Real>& r, std::vector<Real>& z) const {
    // z = M^{-1} * r = r ./ diag(K)
    for (std::size_t i = 0; i < r.size(); ++i) {
        z[i] = diag_inv_[i] * r[i];
    }
}

// ═══════════════════════════════════════════════════════════
// SSOR 预条件器
// ═══════════════════════════════════════════════════════════

void SSORPreconditioner::build(const SparseMatrixCSR& K) {
    K_ = &K;
    std::size_t n = K.rows();
    diag_inv_.resize(n);
    
    // 提取对角元素的倒数
    for (std::size_t i = 0; i < n; ++i) {
        Real diag = 0.0;
        for (std::size_t idx = K.row_ptr()[i]; idx < K.row_ptr()[i + 1]; ++idx) {
            if (K.col_indices()[idx] == i) {
                diag = K.values()[idx];
                break;
            }
        }
        
        if (std::abs(diag) < 1e-30) {
            diag_inv_[i] = 1.0;
        } else {
            diag_inv_[i] = 1.0 / diag;
        }
    }
    
    FEM_INFO("SSOR preconditioner built: " + std::to_string(n) + " DOFs, omega=" + 
             std::to_string(omega_));
}

void SSORPreconditioner::apply(const std::vector<Real>& r, std::vector<Real>& z) const {
    // SSOR 预条件器：M^{-1} * r
    // M = (D + ωL) D^{-1} (D + ωU)
    // 
    // 求解分两步：
    // 1. 前向：(D + ωL) y = r
    // 2. 后向：(D + ωU) z = D*y
    
    std::size_t n = r.size();
    std::vector<Real> y(n);
    
    const auto& row_ptr = K_->row_ptr();
    const auto& col_idx = K_->col_indices();
    const auto& values = K_->values();
    
    // 1. 前向替代：(D + ωL) y = r
    for (std::size_t i = 0; i < n; ++i) {
        Real sum = r[i];
        
        // 减去 L * y (j < i)
        for (std::size_t idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
            std::size_t j = col_idx[idx];
            if (j < i) {
                sum -= omega_ * values[idx] * y[j];
            }
        }
        
        y[i] = sum * diag_inv_[i];
    }
    
    // 2. 后向替代：(D + ωU) z = D*y
    for (std::size_t i = n; i-- > 0;) {
        Real sum = y[i] / diag_inv_[i];  // D * y[i]
        
        // 减去 U * z (j > i)
        for (std::size_t idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
            std::size_t j = col_idx[idx];
            if (j > i) {
                sum -= omega_ * values[idx] * z[j];
            }
        }
        
        z[i] = sum * diag_inv_[i];
    }
}

// ═══════════════════════════════════════════════════════════
// PCG 求解器
// ═══════════════════════════════════════════════════════════

PCGSolver::PCGSolver(const std::string& precond_type, Real omega)
    : precond_type_(precond_type), omega_(omega) {
    
    // 默认参数
    max_iter_ = 10000;
    tol_ = 1e-8;
}

SolveResult PCGSolver::solve(const SparseMatrixCSR& K,
                              const std::vector<Real>& F,
                              std::vector<Real>& x) {
    std::size_t n = F.size();
    x.assign(n, 0.0);
    
    // 构建预条件器
    if (precond_type_ == "jacobi") {
        precond_ = std::make_unique<JacobiPreconditioner>();
    } else if (precond_type_ == "ssor") {
        precond_ = std::make_unique<SSORPreconditioner>(omega_);
    } else {
        // 无预条件器（等价于 CG）
        precond_ = nullptr;
    }
    
    if (precond_) {
        precond_->build(K);
    }
    
    // 初始化
    std::vector<Real> r = F;        // r0 = F - K*x0 = F (因为 x0 = 0)
    std::vector<Real> z(n);         // z = M^{-1} * r
    std::vector<Real> p(n);         // 搜索方向
    std::vector<Real> Ap(n);        // A * p
    
    // z0 = M^{-1} * r0
    if (precond_) {
        precond_->apply(r, z);
    } else {
        z = r;  // 无预条件器
    }
    
    p = z;  // p0 = z0
    
    Real rz = dot(r, z);  // (r, z)
    Real r_norm_init = std::sqrt(dot(r, r));
    
    // 主迭代
    for (std::size_t iter = 0; iter < max_iter_; ++iter) {
        // Ap = K * p
        K.matvec(p.data(), Ap.data());
        
        // α = (r, z) / (p, Ap)
        Real pAp = dot(p, Ap);
        
        if (std::abs(pAp) < 1e-30) {
            FEM_WARN("PCG: (p, Ap) ≈ 0, matrix may not be positive definite");
            return {false, iter, std::sqrt(dot(r, r))};
        }
        
        Real alpha = rz / pAp;
        
        // x_{k+1} = x_k + α * p
        axpy(alpha, p, x);
        
        // r_{k+1} = r_k - α * Ap
        axpy(-alpha, Ap, r);
        
        // 检查收敛
        Real r_norm = std::sqrt(dot(r, r));
        Real rel_res = r_norm / r_norm_init;
        
        if (r_norm < tol_ || rel_res < tol_) {
            FEM_INFO("PCG converged: iter=" + std::to_string(iter + 1) +
                    ", residual=" + fmt_sci(r_norm) +
                    ", rel_res=" + fmt_sci(rel_res));
            return {true, iter + 1, r_norm};
        }
        
        // z_{k+1} = M^{-1} * r_{k+1}
        if (precond_) {
            precond_->apply(r, z);
        } else {
            z = r;
        }
        
        // β = (r_{k+1}, z_{k+1}) / (r_k, z_k)
        Real rz_new = dot(r, z);
        Real beta = rz_new / rz;
        
        // p_{k+1} = z_{k+1} + β * p_k
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        
        rz = rz_new;
    }
    
    FEM_WARN("PCG did not converge in " + std::to_string(max_iter_) + " iterations");
    return {false, max_iter_, std::sqrt(dot(r, r))};
}

}  // namespace fem
