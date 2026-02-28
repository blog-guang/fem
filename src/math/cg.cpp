/**
 * cg.cpp - Conjugate Gradient Solver (Refactored with Vector)
 */

#include "math/cg.h"
#include "math/vector.h"
#include "core/logger.h"

namespace fem {

SolveResult CGSolver::solve(const SparseMatrixCSR& K,
                             const std::vector<Real>& F,
                             std::vector<Real>& x)
{
    std::size_t n = F.size();
    x.assign(n, 0.0);

    // 转换为 Vector
    Vector x_vec(n, 0.0);
    Vector F_vec(F);
    
    // r = F - K*x = F (x0=0)
    Vector r = F_vec;
    Vector p = r;
    Real rr = r.norm_squared();

    Vector Ap(n);
    
    for (std::size_t iter = 0; iter < max_iter_; ++iter) {
        // Ap = K * p
        K.matvec(p.data(), Ap.data());
        
        Real pAp   = p.dot(Ap);
        Real alpha = rr / pAp;

        // 更新解和残差（使用 Vector 运算符）
        x_vec += alpha * p;
        r -= alpha * Ap;

        Real rr_new = r.norm_squared();
        Real res = std::sqrt(rr_new);

        if (res < tol_) {
            // 转换回 std::vector
            x.assign(x_vec.data(), x_vec.data() + n);
            
            FEM_INFO("CG converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(res));
            return {true, iter + 1, res};
        }

        Real beta = rr_new / rr;
        
        // 更新搜索方向（使用 Vector 运算符）
        p = r + beta * p;
        
        rr = rr_new;
    }

    // 转换回 std::vector
    x.assign(x_vec.data(), x_vec.data() + n);
    
    FEM_WARN("CG did not converge in " + std::to_string(max_iter_) + " iterations");
    return {false, max_iter_, std::sqrt(rr)};
}

}  // namespace fem
