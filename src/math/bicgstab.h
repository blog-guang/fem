#pragma once

#include "math/solver.h"
#include "math/vector.h"

namespace fem {

// ── BiCGSTAB (非对称矩阵) ──
class BiCGSTABSolver : public LinearSolver {
public:
    SolveResult solve(const SparseMatrixCSR& K,
                      const Vector& F,
                      Vector& x) override;
};

}  // namespace fem
