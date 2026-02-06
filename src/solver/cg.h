#pragma once

#include "solver/solver.h"

namespace fem {

// ── Conjugate Gradient (对称正定矩阵) ──
class CGSolver : public LinearSolver {
public:
    SolveResult solve(const SparseMatrixCSR& K,
                      const std::vector<Real>& F,
                      std::vector<Real>& x) override;
};

}  // namespace fem
