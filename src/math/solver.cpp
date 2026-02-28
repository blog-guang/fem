#include "math/solver.h"
#include "math/cg.h"
#include "math/bicgstab.h"
#include <stdexcept>

namespace fem {

std::unique_ptr<LinearSolver> create_solver(SolverType type) {
    switch (type) {
        case SolverType::CG:       return std::make_unique<CGSolver>();
        case SolverType::BiCGSTAB: return std::make_unique<BiCGSTABSolver>();
        default:
            throw std::invalid_argument("Unknown solver type");
    }
}

}  // namespace fem
