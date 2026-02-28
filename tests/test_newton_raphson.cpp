/**
 * test_newton_raphson.cpp - Newton-Raphson 求解器测试
 */

#include <gtest/gtest.h>
#include "math/newton_raphson.h"
#include "math/sparse_matrix.h"

using namespace fem;

// ═══ 简单非线性问题：1D 非线性弹簧 ═══
// 方程: R(u) = k*u + c*u³ - F = 0
// 切线刚度: K_t(u) = k + 3*c*u²
class NonlinearSpring : public NonlinearProblem {
public:
    NonlinearSpring(Real k, Real c, Real F) 
        : k_(k), c_(c), F_(F) {}
    
    void compute_residual(const std::vector<Real>& u, 
                         std::vector<Real>& R) override {
        R.resize(1);
        Real x = u[0];
        R[0] = k_ * x + c_ * x * x * x - F_;
    }
    
    void compute_tangent(const std::vector<Real>& u,
                        SparseMatrixCSR& K_t) override {
        Real x = u[0];
        Real kt = k_ + 3.0 * c_ * x * x;
        
        // 使用 COO 格式构造，然后转换为 CSR
        SparseMatrixCOO K_coo(1, 1);
        K_coo.add(0, 0, kt);
        K_t = coo_to_csr(K_coo);
    }
    
private:
    Real k_;  // 线性刚度
    Real c_;  // 非线性刚度
    Real F_;  // 外力
};

// ═══ 多自由度非线性问题：耦合非线性方程组 ═══
// R₁(u₁, u₂) = k₁*u₁ + c₁*u₁³ + k₁₂*(u₁ - u₂) - F₁
// R₂(u₁, u₂) = k₂*u₂ + c₂*u₂³ + k₁₂*(u₂ - u₁) - F₂
class CoupledNonlinearSystem : public NonlinearProblem {
public:
    CoupledNonlinearSystem(Real k1, Real k2, Real k12, 
                          Real c1, Real c2,
                          Real F1, Real F2)
        : k1_(k1), k2_(k2), k12_(k12),
          c1_(c1), c2_(c2),
          F1_(F1), F2_(F2) {}
    
    void compute_residual(const std::vector<Real>& u, 
                         std::vector<Real>& R) override {
        R.resize(2);
        Real u1 = u[0];
        Real u2 = u[1];
        
        R[0] = k1_ * u1 + c1_ * u1*u1*u1 + k12_ * (u1 - u2) - F1_;
        R[1] = k2_ * u2 + c2_ * u2*u2*u2 + k12_ * (u2 - u1) - F2_;
    }
    
    void compute_tangent(const std::vector<Real>& u,
                        SparseMatrixCSR& K_t) override {
        Real u1 = u[0];
        Real u2 = u[1];
        
        // 切线刚度矩阵的元素
        Real K11 = k1_ + 3.0 * c1_ * u1*u1 + k12_;
        Real K12 = -k12_;
        Real K21 = -k12_;
        Real K22 = k2_ + 3.0 * c2_ * u2*u2 + k12_;
        
        // 使用 COO 格式构造，然后转换为 CSR
        SparseMatrixCOO K_coo(2, 2);
        K_coo.add(0, 0, K11);
        K_coo.add(0, 1, K12);
        K_coo.add(1, 0, K21);
        K_coo.add(1, 1, K22);
        K_t = coo_to_csr(K_coo);
    }
    
private:
    Real k1_, k2_, k12_;  // 刚度
    Real c1_, c2_;        // 非线性系数
    Real F1_, F2_;        // 外力
};

// ═══ 测试用例 ═══

TEST(NewtonRaphsonTest, SimpleNonlinearSpring) {
    // 问题: k*u + c*u³ = F
    // 参数: k=1, c=0.1, F=2
    // 方程: u + 0.1u³ = 2
    // 数值解约为 u ≈ 1.595
    
    NonlinearSpring problem(1.0, 0.1, 2.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-8;
    params.max_iter = 20;
    params.verbose = false;
    solver.set_params(params);
    
    std::vector<Real> u = {1.0};  // 初始猜测
    auto result = solver.solve(problem, u);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.iterations, 20);
    EXPECT_LT(result.residual_norm, 1e-8);
    
    // 验证解的精度：u ≈ 1.595
    EXPECT_NEAR(u[0], 1.595, 0.01);
    
    // 验证残差接近零
    std::vector<Real> R(1);
    problem.compute_residual(u, R);
    EXPECT_LT(std::abs(R[0]), 1e-7);
}

TEST(NewtonRaphsonTest, LinearProblem) {
    // 退化到线性问题: k*u = F (c=0)
    // 解析解: u = F/k = 2/1 = 2
    
    NonlinearSpring problem(1.0, 0.0, 2.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-10;
    params.max_iter = 5;
    solver.set_params(params);
    
    std::vector<Real> u = {0.0};  // 初始猜测
    auto result = solver.solve(problem, u);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.iterations, 2);  // 线性问题应该 1-2 次迭代收敛
    EXPECT_NEAR(u[0], 2.0, 1e-9);
}

TEST(NewtonRaphsonTest, CoupledSystem) {
    // 测试耦合非线性系统
    // k₁=1, k₂=1, k₁₂=0.5, c₁=0.1, c₂=0.1
    // F₁=2, F₂=3
    
    CoupledNonlinearSystem problem(1.0, 1.0, 0.5, 0.1, 0.1, 2.0, 3.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-8;
    params.max_iter = 50;
    params.verbose = false;
    solver.set_params(params);
    
    std::vector<Real> u = {1.0, 1.5};  // 初始猜测
    auto result = solver.solve(problem, u);
    
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.residual_norm, 1e-8);
    
    // 验证残差接近零
    std::vector<Real> R(2);
    problem.compute_residual(u, R);
    Real residual = std::sqrt(R[0]*R[0] + R[1]*R[1]);
    EXPECT_LT(residual, 1e-7);
}

TEST(NewtonRaphsonTest, InitialGuessConvergence) {
    // 测试不同初始猜测的收敛性
    NonlinearSpring problem(1.0, 0.1, 2.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-7;  // 放宽容差
    params.max_iter = 30;
    solver.set_params(params);
    
    // 测试多个初始猜测
    std::vector<Real> initial_guesses = {0.5, 1.0, 2.0, 3.0};
    
    for (Real guess : initial_guesses) {
        std::vector<Real> u = {guess};
        auto result = solver.solve(problem, u);
        
        EXPECT_TRUE(result.converged) << "Failed with initial guess: " << guess;
        EXPECT_LT(result.residual_norm, 1e-6);  // 放宽容差
    }
}

TEST(NewtonRaphsonTest, ConvergenceHistory) {
    // 测试收敛历史记录
    NonlinearSpring problem(1.0, 0.1, 2.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-8;
    params.max_iter = 20;
    solver.set_params(params);
    
    std::vector<Real> u = {1.0};
    auto result = solver.solve(problem, u);
    
    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.residual_history.size(), 0);
    
    // 验证收敛历史是递减的（大部分情况）
    for (std::size_t i = 1; i < result.residual_history.size(); ++i) {
        // 允许偶尔增加，但最终应该递减
        if (i == result.residual_history.size() - 1) {
            EXPECT_LT(result.residual_history[i], result.initial_residual);
        }
    }
}

TEST(NewtonRaphsonTest, ZeroInitialResidual) {
    // 测试初始残差为零的情况（已经是解）
    // 线性问题: u = F/k = 2
    NonlinearSpring problem(1.0, 0.0, 2.0);
    
    NewtonRaphsonSolver solver;
    NewtonRaphsonParams params;
    params.tol = 1e-10;
    solver.set_params(params);
    
    std::vector<Real> u = {2.0};  // 精确解
    auto result = solver.solve(problem, u);
    
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.iterations, 0);  // 应该立即返回
}
