#include "solver/bicgstab.h"
#include "core/logger.h"
#include <cmath>

namespace fem {

// 兼容旧代码
using Vector = std::vector<Real>;

static Real dot(const Vector& a, const Vector& b) {
    Real s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

SolveResult BiCGSTABSolver::solve(const CSRMatrix& K,
                                   const std::vector<Real>& F,
                                   std::vector<Real>& x)
{
    std::size_t n = F.size();
    x.assign(n, 0.0);

    Vector r     = F;       // r = b - Ax, x0=0 → r=b
    Vector r_hat = r;       // 任意, 常取 r0
    Real   rho   = 1.0;
    Real   alpha = 1.0;
    Real   omega = 1.0;

    Vector v(n, 0.0), p(n, 0.0), s(n), t(n);

    for (std::size_t iter = 0; iter < max_iter_; ++iter) {
        Real rho_new = dot(r_hat, r);
        if (std::abs(rho_new) < 1e-300) {
            return {false, iter, std::sqrt(dot(r, r))};
        }

        Real beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;

        for (std::size_t i = 0; i < n; ++i)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

        K.matvec(p.data(), v.data());
        alpha = rho / dot(r_hat, v);

        for (std::size_t i = 0; i < n; ++i)
            s[i] = r[i] - alpha * v[i];

        Real s_norm = std::sqrt(dot(s, s));
        if (s_norm < tol_) {
            for (std::size_t i = 0; i < n; ++i) x[i] += alpha * p[i];
            FEM_INFO("BiCGSTAB converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(s_norm));
            return {true, iter + 1, s_norm};
        }

        K.matvec(s.data(), t.data());
        omega = dot(t, s) / dot(t, t);

        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i]  = s[i] - omega * t[i];
        }

        Real res = std::sqrt(dot(r, r));
        if (res < tol_) {
            FEM_INFO("BiCGSTAB converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(res));
            return {true, iter + 1, res};
        }
    }

    FEM_WARN("BiCGSTAB did not converge in " + std::to_string(max_iter_) + " iterations");
    return {false, max_iter_, std::sqrt(dot(r, r))};
}

}  // namespace fem
