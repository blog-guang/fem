#include <gtest/gtest.h>
#include "math/solver.h"
#include "math/sparse_matrix.h"

using namespace fem;

// ── 辅助: 构造 SPD 三对角矩阵 ──
// [[2,-1,0],[-1,2,-1],[0,-1,2]]
static SparseMatrixCSR make_tridiag_3() {
    SparseMatrixCOO coo(3, 3);
    coo.add(0,0, 2); coo.add(0,1,-1);
    coo.add(1,0,-1); coo.add(1,1, 2); coo.add(1,2,-1);
    coo.add(2,1,-1); coo.add(2,2, 2);
    return coo_to_csr(coo);
}

// ── CG ──
TEST(CGTest, Solve_2x2) {
    // [[2,-1],[-1,2]] * x = [1,1] → x = [1,1]
    SparseMatrixCOO coo(0, 0);
    coo = SparseMatrixCOO(2, 2);
    coo.add(0,0, 2); coo.add(0,1,-1);
    coo.add(1,0,-1); coo.add(1,1, 2);
    auto K = coo_to_csr(coo);
    std::vector<Real> F = {1.0, 1.0};
    std::vector<Real> x;

    auto solver = create_solver(SolverType::CG);
    auto res = solver->solve(K, F, x);

    EXPECT_TRUE(res.converged);
    EXPECT_NEAR(x[0], 1.0, 1e-8);
    EXPECT_NEAR(x[1], 1.0, 1e-8);
}

TEST(CGTest, Solve_3x3_Tridiag) {
    auto K = make_tridiag_3();
    std::vector<Real> F = {1.0, 0.0, 1.0};
    std::vector<Real> x;

    auto solver = create_solver(SolverType::CG);
    auto res = solver->solve(K, F, x);

    EXPECT_TRUE(res.converged);
    // 验证 Kx ≈ F
    EXPECT_NEAR(2*x[0] - x[1],           F[0], 1e-8);
    EXPECT_NEAR(-x[0] + 2*x[1] - x[2],  F[1], 1e-8);
    EXPECT_NEAR(-x[1] + 2*x[2],          F[2], 1e-8);
}

TEST(CGTest, Solve_Identity) {
    // I * x = [3, 7] → x = [3, 7]
    SparseMatrixCOO coo_I(2, 2);
    coo_I.add(0, 0, 1.0);
    coo_I.add(1, 1, 1.0);
    SparseMatrixCSR I = coo_to_csr(coo_I);

    std::vector<Real> F = {3.0, 7.0};
    std::vector<Real> x;

    auto solver = create_solver(SolverType::CG);
    auto res = solver->solve(I, F, x);

    EXPECT_TRUE(res.converged);
    EXPECT_NEAR(x[0], 3.0, 1e-8);
    EXPECT_NEAR(x[1], 7.0, 1e-8);
}

// ── BiCGSTAB ──
TEST(BiCGSTABTest, Solve_2x2) {
    SparseMatrixCOO coo(0, 0);
    coo = SparseMatrixCOO(2, 2);
    coo.add(0,0, 2); coo.add(0,1,-1);
    coo.add(1,0,-1); coo.add(1,1, 2);
    auto K = coo_to_csr(coo);
    std::vector<Real> F = {1.0, 1.0};
    std::vector<Real> x;

    auto solver = create_solver(SolverType::BiCGSTAB);
    auto res = solver->solve(K, F, x);

    EXPECT_TRUE(res.converged);
    EXPECT_NEAR(x[0], 1.0, 1e-8);
    EXPECT_NEAR(x[1], 1.0, 1e-8);
}

TEST(BiCGSTABTest, Solve_3x3_NonSymmetric) {
    // 非对称矩阵: [[3,1,0],[0,2,1],[1,0,3]]
    // b = [4, 3, 4] → x = [1, 1, 1]
    SparseMatrixCOO coo(0, 0);
    coo = SparseMatrixCOO(3, 3);
    coo.add(0,0,3); coo.add(0,1,1);
    coo.add(1,1,2); coo.add(1,2,1);
    coo.add(2,0,1); coo.add(2,2,3);
    auto K = coo_to_csr(coo);
    std::vector<Real> F = {4.0, 3.0, 4.0};
    std::vector<Real> x;

    auto solver = create_solver(SolverType::BiCGSTAB);
    auto res = solver->solve(K, F, x);

    EXPECT_TRUE(res.converged);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
    EXPECT_NEAR(x[1], 1.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);
}
