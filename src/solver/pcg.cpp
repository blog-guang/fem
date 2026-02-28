/**
 * pcg.cpp - PCG Solver Implementation
 */

#include "solver/pcg.h"
#include "core/logger.h"
#include <cmath>
#include <algorithm>

// AMGCL headers (仅在 .cpp 中包含)
// 定义 AMGCL_NO_BOOST 以避免 Boost 依赖
#define AMGCL_NO_BOOST
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

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
// ILU(0) 预条件器
// ═══════════════════════════════════════════════════════════

/**
 * 在 CSR 矩阵的第 row 行中查找列 col 的索引
 * 返回 idx 使得 col_indices[idx] == col
 * 如果不存在，返回 -1
 */
static Index find_col_in_row(const SparseMatrixCSR& mat, Index row, Index col) {
    const auto& row_ptr = mat.row_ptr();
    const auto& col_idx = mat.col_indices();
    
    for (Index idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx) {
        if (col_idx[idx] == col) {
            return idx;
        }
    }
    return -1;  // 未找到
}

void ILUPreconditioner::build(const SparseMatrixCSR& K) {
    // 1. 复制 K 到 LU_（将原位修改为 L*U）
    LU_ = K;
    
    std::size_t n = LU_.rows();
    auto& values = const_cast<std::vector<Real>&>(LU_.values());
    const auto& row_ptr = LU_.row_ptr();
    const auto& col_idx = LU_.col_indices();
    
    // 2. ILU(0) 分解
    // 对每一行 i 进行消元
    for (std::size_t i = 0; i < n; ++i) {
        // 找到第 i 行的对角元素
        Index diag_idx = find_col_in_row(LU_, i, i);
        
        if (diag_idx < 0 || std::abs(values[diag_idx]) < 1e-30) {
            FEM_WARN("ILU(0): row " + std::to_string(i) + 
                    " has near-zero diagonal, using 1.0");
            if (diag_idx >= 0) {
                values[diag_idx] = 1.0;
            }
            continue;
        }
        
        Real diag = values[diag_idx];
        
        // 对第 i 行中所有列 k < i 的元素（下三角部分）
        for (Index k_idx = row_ptr[i]; k_idx < row_ptr[i + 1]; ++k_idx) {
            Index k = col_idx[k_idx];
            
            if (k >= i) break;  // 只处理下三角
            
            // L(i, k) = A(i, k) / U(k, k)
            Index k_diag_idx = find_col_in_row(LU_, k, k);
            
            if (k_diag_idx < 0 || std::abs(values[k_diag_idx]) < 1e-30) {
                continue;  // 跳过病态行
            }
            
            values[k_idx] /= values[k_diag_idx];
            Real L_ik = values[k_idx];
            
            // 更新 U 部分：A(i, j) -= L(i, k) * U(k, j)
            // 对第 i 行中所有列 j > k 的元素
            for (Index j_idx = row_ptr[i]; j_idx < row_ptr[i + 1]; ++j_idx) {
                Index j = col_idx[j_idx];
                
                if (j <= k) continue;  // 只更新 j > k 的部分
                
                // 查找 U(k, j)
                Index kj_idx = find_col_in_row(LU_, k, j);
                
                if (kj_idx >= 0) {
                    // A(i, j) -= L(i, k) * U(k, j)
                    values[j_idx] -= L_ik * values[kj_idx];
                }
            }
        }
    }
    
    FEM_INFO("ILU(0) preconditioner built: " + std::to_string(n) + " DOFs");
}

void ILUPreconditioner::forward_solve(const std::vector<Real>& r, 
                                       std::vector<Real>& y) const {
    // 求解 L*y = r
    // L 是单位下三角矩阵（对角线为 1）
    
    std::size_t n = r.size();
    y.resize(n);
    
    const auto& row_ptr = LU_.row_ptr();
    const auto& col_idx = LU_.col_indices();
    const auto& values = LU_.values();
    
    for (std::size_t i = 0; i < n; ++i) {
        Real sum = r[i];
        
        // sum -= L(i, j) * y[j]  (j < i)
        for (Index idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
            Index j = col_idx[idx];
            
            if (j >= i) break;  // 只处理下三角
            
            sum -= values[idx] * y[j];
        }
        
        y[i] = sum;  // L 的对角线为 1
    }
}

void ILUPreconditioner::backward_solve(const std::vector<Real>& y, 
                                        std::vector<Real>& z) const {
    // 求解 U*z = y
    // U 是上三角矩阵
    
    std::size_t n = y.size();
    z.resize(n);
    
    const auto& row_ptr = LU_.row_ptr();
    const auto& col_idx = LU_.col_indices();
    const auto& values = LU_.values();
    
    for (std::size_t i = n; i-- > 0;) {
        Real sum = y[i];
        Real diag = 1.0;
        
        // sum -= U(i, j) * z[j]  (j > i)
        // 同时找到对角元素 U(i, i)
        for (Index idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
            Index j = col_idx[idx];
            
            if (j == i) {
                diag = values[idx];
            } else if (j > i) {
                sum -= values[idx] * z[j];
            }
        }
        
        if (std::abs(diag) < 1e-30) {
            z[i] = sum;  // 避免除零
        } else {
            z[i] = sum / diag;
        }
    }
}

void ILUPreconditioner::apply(const std::vector<Real>& r, std::vector<Real>& z) const {
    // 求解 M^{-1} * r，其中 M = L * U
    // 
    // 分两步：
    // 1. 前向替换: L * y = r
    // 2. 后向替换: U * z = y
    
    std::vector<Real> y;
    forward_solve(r, y);
    backward_solve(y, z);
}

// ═══════════════════════════════════════════════════════════
// AMG 预条件器 (AMGCL)
// ═══════════════════════════════════════════════════════════

/**
 * AMG 预条件器的内部实现类 (PIMPL)
 * 
 * 使用 AMGCL 的内置后端和 smoothed aggregation 粗化策略
 */
class AMGPreconditioner::Impl {
public:
    // AMGCL 类型定义
    typedef amgcl::backend::builtin<double> Backend;
    
    typedef amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    > AMG;
    
    std::unique_ptr<AMG> amg;
    std::size_t n;
};

AMGPreconditioner::AMGPreconditioner()
    : pimpl_(std::make_unique<Impl>()) {
}

AMGPreconditioner::~AMGPreconditioner() = default;

void AMGPreconditioner::build(const SparseMatrixCSR& K) {
    pimpl_->n = K.rows();
    
    // 将 CSR 矩阵转换为 AMGCL 格式
    // AMGCL 使用 boost::tie(n, ptr, col, val) 作为输入
    // 但不使用 Boost 时，直接使用 adapter
    
    // 注意：AMGCL 需要 ptrdiff_t 类型的指针数组
    std::vector<ptrdiff_t> ptr(K.row_ptr().begin(), K.row_ptr().end());
    std::vector<ptrdiff_t> col(K.col_indices().begin(), K.col_indices().end());
    std::vector<double> val(K.values().begin(), K.values().end());
    
    // 构建 AMG 层次结构
    // 不使用 Boost property_tree，直接设置参数结构体
    typename Impl::AMG::params prm;
    prm.coarsening.aggr.eps_strong = 0.0;  // 强连接阈值
    
    // 使用 CRS tuple adapter
    pimpl_->amg = std::make_unique<Impl::AMG>(
        std::tie(pimpl_->n, ptr, col, val), 
        prm
    );
    
    FEM_INFO("AMG preconditioner built: " + std::to_string(pimpl_->n) + " DOFs");
}

void AMGPreconditioner::apply(const std::vector<Real>& r, std::vector<Real>& z) const {
    // 使用 AMG V-cycle 作为预条件器
    // z = M^{-1} * r
    
    z.resize(r.size());
    std::fill(z.begin(), z.end(), 0.0);
    
    // AMGCL 的 apply 方法：执行一次 V-cycle
    pimpl_->amg->apply(r, z);
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
    } else if (precond_type_ == "ilu") {
        precond_ = std::make_unique<ILUPreconditioner>();
    } else if (precond_type_ == "amg") {
        precond_ = std::make_unique<AMGPreconditioner>();
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
