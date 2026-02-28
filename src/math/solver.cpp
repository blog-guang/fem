#include "math/solver.h"
#include "math/cg.h"
#include "math/bicgstab.h"
#include <stdexcept>

namespace fem {

// 兼容旧接口：std::vector → Vector 转换
SolveResult LinearSolver::solve(const SparseMatrixCSR& K,
                                 const std::vector<Real>& F,
                                 std::vector<Real>& x) {
    // 转换为 Vector
    Vector F_vec(F);
    Vector x_vec(x.size());
    
    // 调用新接口
    SolveResult result = solve(K, F_vec, x_vec);
    
    // 转换回 std::vector
    x.assign(x_vec.data(), x_vec.data() + x_vec.size());
    
    return result;
}

std::unique_ptr<LinearSolver> create_solver(SolverType type) {
    switch (type) {
        case SolverType::CG:       return std::make_unique<CGSolver>();
        case SolverType::BiCGSTAB: return std::make_unique<BiCGSTABSolver>();
        default:
            throw std::invalid_argument("Unknown solver type");
    }
}

}  // namespace fem
