#pragma once

#include "core/types.h"
#include "math/sparse_matrix.h"
#include <memory>

namespace fem {

// ── 求解结果 ──
struct SolveResult {
    bool        converged;
    std::size_t iterations;
    Real        residual;      // 最终残差 ||r||
};

// ── 线性求解器基类 ──
// 求解 Kx = F
class LinearSolver {
public:
    virtual ~LinearSolver() = default;

    virtual SolveResult solve(const SparseMatrixCSR& K,
                              const std::vector<Real>& F,
                              std::vector<Real>& x) = 0;

    void set_tol(Real tol)              { tol_ = tol; }
    void set_max_iter(std::size_t iter) { max_iter_ = iter; }

protected:
    Real        tol_     {1e-10};
    std::size_t max_iter_{1000};
};

// ── 工厂 ──
enum class SolverType { CG, BiCGSTAB };

[[nodiscard]] std::unique_ptr<LinearSolver> create_solver(SolverType type);

}  // namespace fem
