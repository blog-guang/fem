#pragma once

#include "solver/solver.h"

namespace fem {

// ── BiCGSTAB (非对称矩阵) ──
class BiCGSTABSolver : public LinearSolver {
public:
    SolveResult solve(const CSRMatrix& K,
                      const Vector&    F,
                      Vector&          x) override;
};

}  // namespace fem
