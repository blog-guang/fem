/**
 * test_amg.cpp - AMG Preconditioner Tests (using AMGCL)
 */

#include <gtest/gtest.h>
#include "math/pcg.h"
#include "math/sparse_matrix.h"
#include <vector>

using namespace fem;

/**
 * 测试 1: 简单 5x5 Poisson 问题
 */
TEST(AMGTest, Simple5x5) {
    // 5 节点 1D Poisson 问题
    // K = [ 2  -1   0   0   0 ]
    //     [-1   2  -1   0   0 ]
    //     [ 0  -1   2  -1   0 ]
    //     [ 0   0  -1   2  -1 ]
    //     [ 0   0   0  -1   2 ]
    
    SparseMatrixCSR K;
    std::vector<Index> row_ptr = {0, 2, 5, 8, 11, 13};
    std::vector<Index> col_idx = {0, 1,  0, 1, 2,  1, 2, 3,  2, 3, 4,  3, 4};
    std::vector<Real> values   = {2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1, 2};
    
    K.set_data(5, 5, std::move(row_ptr), std::move(col_idx), std::move(values));
    
    std::vector<Real> F = {1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<Real> x(5);
    
    // 使用 AMG 预条件器求解
    PCGSolver solver_amg("amg");
    solver_amg.set_tol(1e-8);
    
    auto result = solver_amg.solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual, 1e-8);
    
    // 验证解的合理性
    for (Real xi : x) {
        EXPECT_FALSE(std::isnan(xi));
        EXPECT_FALSE(std::isinf(xi));
    }
}

/**
 * 测试 2: AMG vs ILU vs Jacobi 收敛速度对比
 */
TEST(AMGTest, CompareWithOthers) {
    // 10x10 三对角矩阵
    std::size_t n = 10;
    SparseMatrixCOO coo(n, n);
    
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) coo.add(i, i - 1, -1.0);
        coo.add(i, i, 2.0);
        if (i < n - 1) coo.add(i, i + 1, -1.0);
    }
    
    SparseMatrixCSR K = coo_to_csr(coo);
    std::vector<Real> F(n, 1.0);
    
    std::vector<Real> x_jacobi(n), x_ilu(n), x_amg(n);
    
    // Jacobi
    PCGSolver solver_jacobi("jacobi");
    solver_jacobi.set_tol(1e-8);
    auto result_jacobi = solver_jacobi.solve(K, F, x_jacobi);
    
    // ILU
    PCGSolver solver_ilu("ilu");
    solver_ilu.set_tol(1e-8);
    auto result_ilu = solver_ilu.solve(K, F, x_ilu);
    
    // AMG
    PCGSolver solver_amg("amg");
    solver_amg.set_tol(1e-8);
    auto result_amg = solver_amg.solve(K, F, x_amg);
    
    // 都应该收敛
    EXPECT_TRUE(result_jacobi.converged);
    EXPECT_TRUE(result_ilu.converged);
    EXPECT_TRUE(result_amg.converged);
    
    // AMG 应该很快收敛
    EXPECT_LT(result_amg.iterations, result_jacobi.iterations);
    
    std::cout << "Jacobi: " << result_jacobi.iterations << " iterations\n";
    std::cout << "ILU(0): " << result_ilu.iterations << " iterations\n";
    std::cout << "AMG:    " << result_amg.iterations << " iterations\n";
}

/**
 * 测试 3: 大型稀疏矩阵 (50x50)
 */
TEST(AMGTest, LargeSystem) {
    std::size_t n = 50;
    SparseMatrixCOO coo(n, n);
    
    // 三对角矩阵
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) coo.add(i, i - 1, -1.0);
        coo.add(i, i, 3.0);  // 对角占优
        if (i < n - 1) coo.add(i, i + 1, -1.0);
    }
    
    SparseMatrixCSR K = coo_to_csr(coo);
    std::vector<Real> F(n, 1.0);
    std::vector<Real> x(n);
    
    PCGSolver solver("amg");
    solver.set_tol(1e-10);
    solver.set_max_iter(1000);
    
    auto result = solver.solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual, 1e-10);
    EXPECT_LT(result.iterations, 50);  // AMG 应该很快收敛
}

/**
 * 测试 4: 2D Poisson 问题（更接近实际FEM问题）
 */
TEST(AMGTest, Poisson2D) {
    // 生成 10x10 2D Poisson 问题的五点差分矩阵
    int nx = 10, ny = 10;
    int n = nx * ny;
    
    SparseMatrixCOO coo(n, n);
    
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            int idx = i * nx + j;
            
            // 对角元素
            coo.add(idx, idx, 4.0);
            
            // 左
            if (j > 0) coo.add(idx, idx - 1, -1.0);
            
            // 右
            if (j < nx - 1) coo.add(idx, idx + 1, -1.0);
            
            // 下
            if (i > 0) coo.add(idx, idx - nx, -1.0);
            
            // 上
            if (i < ny - 1) coo.add(idx, idx + nx, -1.0);
        }
    }
    
    SparseMatrixCSR K = coo_to_csr(coo);
    std::vector<Real> F(n, 1.0);
    std::vector<Real> x(n);
    
    PCGSolver solver("amg");
    solver.set_tol(1e-8);
    solver.set_max_iter(1000);
    
    auto result = solver.solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual, 1e-8);
    EXPECT_LT(result.iterations, 30);  // AMG 对 2D Poisson 应该非常高效
    
    std::cout << "2D Poisson (10x10): " << result.iterations << " iterations\n";
}
