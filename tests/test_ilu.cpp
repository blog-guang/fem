/**
 * test_ilu.cpp - ILU(0) Preconditioner Tests
 */

#include <gtest/gtest.h>
#include "math/pcg.h"
#include "math/sparse_matrix.h"
#include <vector>

using namespace fem;

/**
 * 测试 1: 简单 3x3 对称正定矩阵
 * 
 * K = [ 4  -1   0 ]
 *     [-1   4  -1 ]
 *     [ 0  -1   4 ]
 * 
 * 验证 ILU(0) 分解的正确性
 */
TEST(ILUTest, SimpleMatrix3x3) {
    // 构造 CSR 矩阵
    SparseMatrixCSR K;
    std::vector<Index> row_ptr = {0, 2, 5, 7};
    std::vector<Index> col_idx = {0, 1,  0, 1, 2,  1, 2};
    std::vector<Real> values   = {4, -1, -1, 4, -1, -1, 4};
    
    K.set_data(3, 3, std::move(row_ptr), std::move(col_idx), std::move(values));
    
    // 构建 ILU(0) 预条件器
    ILUPreconditioner ilu;
    ilu.build(K);
    
    // 测试预条件器应用: z = M^{-1} * r
    std::vector<Real> r = {1.0, 2.0, 3.0};
    std::vector<Real> z(3);
    
    ilu.apply(r, z);
    
    // z 应该接近 K^{-1} * r（不完全相等，因为是不完全分解）
    // 至少验证 z 不是零向量
    EXPECT_GT(std::abs(z[0]), 1e-10);
    EXPECT_GT(std::abs(z[1]), 1e-10);
    EXPECT_GT(std::abs(z[2]), 1e-10);
}

/**
 * 测试 2: 对角占优矩阵（ILU 稳定的情况）
 */
TEST(ILUTest, DiagonallyDominant) {
    // K = [ 10  -1   0 ]
    //     [ -1  10  -1 ]
    //     [  0  -1  10 ]
    
    SparseMatrixCSR K;
    std::vector<Index> row_ptr = {0, 2, 5, 7};
    std::vector<Index> col_idx = {0, 1,  0, 1, 2,  1, 2};
    std::vector<Real> values   = {10, -1, -1, 10, -1, -1, 10};
    
    K.set_data(3, 3, std::move(row_ptr), std::move(col_idx), std::move(values));
    
    ILUPreconditioner ilu;
    ilu.build(K);
    
    std::vector<Real> r = {1.0, 0.0, 0.0};
    std::vector<Real> z(3);
    
    ilu.apply(r, z);
    
    // 对角占优矩阵的 ILU(0) 应该很稳定
    EXPECT_FALSE(std::isnan(z[0]));
    EXPECT_FALSE(std::isnan(z[1]));
    EXPECT_FALSE(std::isnan(z[2]));
}

/**
 * 测试 3: PCG 求解器 + ILU 预条件器
 * 
 * 求解 Poisson 方程的离散系统
 */
TEST(ILUTest, PCGSolveWithILU) {
    // 5 节点 1D Poisson 问题的刚度矩阵
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
    Vector x(5);
    
    // 使用 ILU 预条件器求解
    PCGSolver solver_ilu("ilu");
    solver_ilu.set_tol(1e-8);
    
    auto result = solver_ilu.solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual, 1e-8);
    
    // 验证解的合理性
    for (Real xi : x) {
        EXPECT_FALSE(std::isnan(xi));
        EXPECT_FALSE(std::isinf(xi));
    }
}

/**
 * 测试 4: ILU vs Jacobi 收敛速度对比
 */
TEST(ILUTest, CompareWithJacobi) {
    // 较大的系统 (10x10)
    std::size_t n = 10;
    std::vector<Index> row_ptr(n + 1);
    std::vector<Index> col_idx;
    std::vector<Real> values;
    
    row_ptr[0] = 0;
    
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) {
            col_idx.push_back(i - 1);
            values.push_back(-1.0);
        }
        
        col_idx.push_back(i);
        values.push_back(2.0);
        
        if (i < n - 1) {
            col_idx.push_back(i + 1);
            values.push_back(-1.0);
        }
        
        row_ptr[i + 1] = col_idx.size();
    }
    
    SparseMatrixCSR K;
    K.set_data(n, n, std::move(row_ptr), std::move(col_idx), std::move(values));
    
    std::vector<Real> F(n, 1.0);
    Vector x_jacobi(n), x_ilu(n);
    
    // Jacobi 预条件
    PCGSolver solver_jacobi("jacobi");
    solver_jacobi.set_tol(1e-8);
    auto result_jacobi = solver_jacobi.solve(K, F, x_jacobi);
    
    // ILU 预条件
    PCGSolver solver_ilu("ilu");
    solver_ilu.set_tol(1e-8);
    auto result_ilu = solver_ilu.solve(K, F, x_ilu);
    
    // 两者都应该收敛
    EXPECT_TRUE(result_jacobi.converged);
    EXPECT_TRUE(result_ilu.converged);
    
    // ILU 应该更快收敛（迭代次数更少）
    EXPECT_LT(result_ilu.iterations, result_jacobi.iterations);
    
    std::cout << "Jacobi: " << result_jacobi.iterations << " iterations\n";
    std::cout << "ILU(0): " << result_ilu.iterations << " iterations\n";
    std::cout << "Speedup: " << static_cast<double>(result_jacobi.iterations) / result_ilu.iterations << "x\n";
}

/**
 * 测试 5: 大型稀疏矩阵 (50x50)
 */
TEST(ILUTest, LargeSystem) {
    std::size_t n = 50;
    std::vector<Index> row_ptr(n + 1);
    std::vector<Index> col_idx;
    std::vector<Real> values;
    
    row_ptr[0] = 0;
    
    // 三对角矩阵
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) {
            col_idx.push_back(i - 1);
            values.push_back(-1.0);
        }
        
        col_idx.push_back(i);
        values.push_back(3.0);  // 对角占优
        
        if (i < n - 1) {
            col_idx.push_back(i + 1);
            values.push_back(-1.0);
        }
        
        row_ptr[i + 1] = col_idx.size();
    }
    
    SparseMatrixCSR K;
    K.set_data(n, n, std::move(row_ptr), std::move(col_idx), std::move(values));
    
    std::vector<Real> F(n, 1.0);
    Vector x(n);
    
    PCGSolver solver("ilu");
    solver.set_tol(1e-10);
    solver.set_max_iter(1000);
    
    auto result = solver.solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual, 1e-10);
    EXPECT_LT(result.iterations, 100);  // ILU 应该很快收敛
}
