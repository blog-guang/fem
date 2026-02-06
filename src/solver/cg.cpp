#include "solver/cg.h"
#include "core/logger.h"
#include <cmath>

namespace fem {

static Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

SolveResult CGSolver::solve(const SparseMatrixCSR& K,
                             const std::vector<Real>& F,
                             std::vector<Real>& x)
{
    std::size_t n = F.size();
    x.assign(n, 0.0);

    // r = F - K*x = F (x0=0)
    std::vector<Real> r = F;
    std::vector<Real> p = r;
    Real   rr = dot(r, r);

    std::vector<Real> Ap(n);

    for (std::size_t iter = 0; iter < max_iter_; ++iter) {
        K.matvec(p.data(), Ap.data());
        Real pAp   = dot(p, Ap);
        Real alpha = rr / pAp;

        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        Real rr_new = dot(r, r);
        Real res    = std::sqrt(rr_new);

        if (res < tol_) {
            FEM_INFO("CG converged: iter=" + std::to_string(iter + 1) +
                     " residual=" + fmt_sci(res));
            return {true, iter + 1, res};
        }

        Real beta = rr_new / rr;
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rr = rr_new;
    }

    FEM_WARN("CG did not converge in " + std::to_string(max_iter_) + " iterations");
    return {false, max_iter_, std::sqrt(rr)};
}

}  // namespace fem
