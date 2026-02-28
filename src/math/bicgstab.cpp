/**
 * bicgstab.cpp - BiCGSTAB Solver (Refactored with Vector)
 */

#include "math/bicgstab.h"
#include "math/vector.h"
#include "core/logger.h"

namespace fem {

SolveResult BiCGSTABSolver::solve(const SparseMatrixCSR& K,
                                   const Vector& F,
                                   Vector& x)
{
    std::size_t n = F.size();
    x = Vector(n, 0.0);

    Vector r     = F;               // r = b - Ax, x0=0 → r=b
    Vector r_hat = r;               // 任意, 常取 r0
    Real   rho   = 1.0;
    Real   alpha = 1.0;
    Real   omega = 1.0;

    Vector v(n, 0.0), p(n, 0.0), s(n), t(n);

    for (std::size_t iter = 0; iter < max_iter_; ++iter) {
        Real rho_new = r_hat.dot(r);
        if (std::abs(rho_new) < 1e-300) {
            return {false, iter, r.norm()};
        }

        Real beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;

        // p = r + β * (p - ω * v)
        p = r + beta * (p - omega * v);

        // v = K * p
        K.matvec(p.data(), v.data());
        
        alpha = rho / r_hat.dot(v);

        // s = r - α * v
        s = r - alpha * v;

        Real s_norm = s.norm();
        if (s_norm < tol_) {
            x += alpha * p;
            FEM_INFO("BiCGSTAB converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(s_norm));
            return {true, iter + 1, s_norm};
        }

        // t = K * s
        K.matvec(s.data(), t.data());
        
        omega = t.dot(s) / t.dot(t);

        // x += α*p + ω*s
        x += alpha * p + omega * s;
        
        // r = s - ω*t
        r = s - omega * t;

        Real res = r.norm();
        if (res < tol_) {
            FEM_INFO("BiCGSTAB converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(res));
            return {true, iter + 1, res};
        }
    }

    FEM_WARN("BiCGSTAB did not converge in " + std::to_string(max_iter_) + " iterations");
    return {false, max_iter_, r.norm()};
}

}  // namespace fem
