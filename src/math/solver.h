#pragma once

#include "core/types.h"
#include "math/sparse_matrix.h"
#include "math/vector.h"
#include <memory>
#include <vector>

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

    // 新接口：使用 Vector
    virtual SolveResult solve(const SparseMatrixCSR& K,
                              const Vector& F,
                              Vector& x) = 0;
    
    // 兼容旧接口：std::vector（内部转换为 Vector）
    SolveResult solve(const SparseMatrixCSR& K,
                      const std::vector<Real>& F,
                      std::vector<Real>& x);

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
