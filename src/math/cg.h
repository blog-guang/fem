#pragma once

#include "math/solver.h"
#include "math/vector.h"

namespace fem {

// ── Conjugate Gradient (对称正定矩阵) ──
class CGSolver : public LinearSolver {
public:
    SolveResult solve(const SparseMatrixCSR& K,
                      const Vector& F,
                      Vector& x) override;
};

}  // namespace fem
