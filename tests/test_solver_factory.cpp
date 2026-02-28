#include <gtest/gtest.h>
#include "math/solver_factory.h"
#include "math/cg.h"
#include "math/sparse_matrix.h"
#include "math/vector.h"

using namespace fem;

// ═══════════════════════════════════════════════════════════
// SolverFactory 测试
// ═══════════════════════════════════════════════════════════

TEST(SolverFactoryTest, CreateCG) {
    auto solver = SolverFactory::create("CG", {
        {"tol", 1e-8},
        {"max_iter", 500}
    });
    
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, CreatePCG) {
    auto solver = SolverFactory::create("PCG", 
        {{"tol", 1e-10}},
        {{"precond", "jacobi"}}
    );
    
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, CreatePCGWithSSOR) {
    auto solver = SolverFactory::create("PCG",
        {{"tol", 1e-10}, {"omega", 1.2}},
        {{"precond", "ssor"}}
    );
    
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, CreateBiCGSTAB) {
    auto solver = SolverFactory::create("BiCGSTAB", {
        {"tol", 1e-8}
    });
    
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, UnknownType) {
    EXPECT_THROW(
        SolverFactory::create("NonExistentSolver"),
        std::invalid_argument
    );
}

TEST(SolverFactoryTest, IsRegistered) {
    EXPECT_TRUE(SolverFactory::isRegistered("CG"));
    EXPECT_TRUE(SolverFactory::isRegistered("PCG"));
    EXPECT_TRUE(SolverFactory::isRegistered("BiCGSTAB"));
    EXPECT_FALSE(SolverFactory::isRegistered("UnknownSolver"));
}

TEST(SolverFactoryTest, GetRegisteredTypes) {
    auto types = SolverFactory::getRegisteredTypes();
    
    EXPECT_GE(types.size(), 3);  // 至少 CG, PCG, BiCGSTAB
    
    bool found_cg = false;
    bool found_pcg = false;
    bool found_bicgstab = false;
    
    for (const auto& name : types) {
        if (name == "CG") found_cg = true;
        if (name == "PCG") found_pcg = true;
        if (name == "BiCGSTAB") found_bicgstab = true;
    }
    
    EXPECT_TRUE(found_cg);
    EXPECT_TRUE(found_pcg);
    EXPECT_TRUE(found_bicgstab);
}

TEST(SolverFactoryTest, DefaultParameters) {
    // 使用默认参数创建（空参数列表）
    auto solver = SolverFactory::create("CG");
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, RegisterCustomSolver) {
    // 注册一个自定义求解器（实际上创建 CG）
    SolverFactory::registerCreator("TestSolver",
        [](const SolverFactory::Parameters& num_p,
           const SolverFactory::StringParams& str_p) 
           -> std::unique_ptr<LinearSolver> 
        {
            return std::make_unique<CGSolver>();
        }
    );
    
    EXPECT_TRUE(SolverFactory::isRegistered("TestSolver"));
    
    auto solver = SolverFactory::create("TestSolver");
    ASSERT_NE(solver, nullptr);
}

TEST(SolverFactoryTest, SolveSimpleSystem) {
    // 通过工厂创建求解器，测试求解功能
    auto solver = SolverFactory::create("CG", {
        {"tol", 1e-8},
        {"max_iter", 100}
    });
    
    // 创建简单的 2x2 系统：K = [[2, 1], [1, 2]]
    SparseMatrixCOO coo(2, 2);
    coo.add(0, 0, 2.0);
    coo.add(0, 1, 1.0);
    coo.add(1, 0, 1.0);
    coo.add(1, 1, 2.0);
    SparseMatrixCSR K = coo_to_csr(coo);
    
    Vector F(2);
    F[0] = 3.0;  // 2*x + y = 3
    F[1] = 3.0;  // x + 2*y = 3
    
    Vector x(2, 0.0);
    
    auto result = solver->solve(K, F, x);
    
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(x[0], 1.0, 1e-6);  // x = 1
    EXPECT_NEAR(x[1], 1.0, 1e-6);  // y = 1
}

TEST(SolverFactoryTest, PCGWithDifferentPreconditioners) {
    // 创建测试矩阵（5x5 对角占优）
    SparseMatrixCOO coo(5, 5);
    for (std::size_t i = 0; i < 5; ++i) {
        coo.add(i, i, 10.0);
        if (i > 0) coo.add(i, i-1, -1.0);
        if (i < 4) coo.add(i, i+1, -1.0);
    }
    SparseMatrixCSR K = coo_to_csr(coo);
    
    Vector F(5, 1.0);
    
    // 测试不同的预条件器
    std::vector<std::string> precond_types = {"jacobi", "ssor", "ilu"};
    
    for (const auto& precond : precond_types) {
        auto solver = SolverFactory::create("PCG",
            {{"tol", 1e-8}},
            {{"precond", precond}}
        );
        
        Vector x(5, 0.0);
        auto result = solver->solve(K, F, x);
        
        EXPECT_TRUE(result.converged) << "Failed with precond: " << precond;
    }
}
