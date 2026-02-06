#pragma once

#include "solver/solver.h"

namespace fem {

// ── BiCGSTAB (非对称矩阵) ──
class BiCGSTABSolver : public LinearSolver {
public:
    SolveResult solve(const CSRMatrix& K,
                      const std::vector<Real>& F,
                      std::vector<Real>& x) override;
};

}  // namespace fem
