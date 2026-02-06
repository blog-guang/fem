#pragma once

#include "solver/solver.h"

namespace fem {

// ── Conjugate Gradient (对称正定矩阵) ──
class CGSolver : public LinearSolver {
public:
    SolveResult solve(const CSRMatrix& K,
                      const Vector&    F,
                      Vector&          x) override;
};

}  // namespace fem
