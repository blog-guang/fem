#include <gtest/gtest.h>
#include "assembly/sparse_matrix.h"
#include "assembly/assembler.h"
#include "assembly/boundary_condition.h"
#include "mesh/mesh_generator.h"

using namespace fem;

// ── SparseMatrix ──
TEST(COOTest, AddAndSize) {
    COOMatrix K;
    K.rows = 3; K.cols = 3;
    K.add(0, 0, 1.0);
    K.add(0, 1, 2.0);
    K.add(1, 1, 3.0);

    EXPECT_EQ(K.row_idx.size(), 3);
    EXPECT_EQ(K.values[2], 3.0);
}

TEST(CSRTest, Matvec_Identity) {
    // 3x3 单位矩阵
    CSRMatrix I;
    I.rows    = 3;
    I.row_ptr = {0, 1, 2, 3};
    I.col_idx = {0, 1, 2};
    I.values  = {1.0, 1.0, 1.0};

    Real x[] = {2.0, 3.0, 5.0};
    Real y[3];
    I.matvec(x, y);

    EXPECT_DOUBLE_EQ(y[0], 2.0);
    EXPECT_DOUBLE_EQ(y[1], 3.0);
    EXPECT_DOUBLE_EQ(y[2], 5.0);
}

TEST(CSRTest, COO_to_CSR_basic) {
    // 2x2: [[2, -1], [-1, 2]]
    COOMatrix coo;
    coo.rows = 2; coo.cols = 2;
    coo.add(0, 0,  2.0);
    coo.add(0, 1, -1.0);
    coo.add(1, 0, -1.0);
    coo.add(1, 1,  2.0);

    auto csr = coo_to_csr(coo);

    EXPECT_EQ(csr.rows, 2);
    // matvec: [2,-1;-1,2] * [1,1] = [1,1]
    Real x[] = {1.0, 1.0};
    Real y[2];
    csr.matvec(x, y);
    EXPECT_NEAR(y[0], 1.0, 1e-15);
    EXPECT_NEAR(y[1], 1.0, 1e-15);
}

TEST(CSRTest, COO_to_CSR_duplicate_sum) {
    // 重复项求和: (0,0) 添加两次 1.0 → 结果 2.0
    COOMatrix coo;
    coo.rows = 2; coo.cols = 2;
    coo.add(0, 0, 1.0);
    coo.add(0, 0, 1.0);  // 重复
    coo.add(1, 1, 3.0);

    auto csr = coo_to_csr(coo);

    Real x[] = {1.0, 1.0};
    Real y[2];
    csr.matvec(x, y);
    EXPECT_NEAR(y[0], 2.0, 1e-15);  // 1+1 = 2
    EXPECT_NEAR(y[1], 3.0, 1e-15);
}

// ── Assembler (简单常数刚度矩阵验证) ──
// 用一个 1x1 三角形网格, 单元刚度 = 单位矩阵 * 面积
static void unit_stiffness(const Vec3* /*coords*/, std::size_t n_nodes, std::size_t dofs_per_node, Real* Ke, void* /*ctx*/) {
    std::size_t n = n_nodes * dofs_per_node;
    for (std::size_t i = 0; i < n * n; ++i) Ke[i] = 0.0;
    for (std::size_t i = 0; i < n; ++i)     Ke[i * n + i] = 1.0;
}

static void unit_load(const Vec3* /*coords*/, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* /*ctx*/) {
    std::size_t n = n_nodes * dofs_per_node;
    for (std::size_t i = 0; i < n; ++i) Fe[i] = 1.0;
}

TEST(AssemblerTest, SingleTriangle) {
    // 手动 1 个三角形
    Mesh mesh;
    mesh.add_node({0.0, 0.0, 0.0});
    mesh.add_node({1.0, 0.0, 0.0});
    mesh.add_node({0.0, 1.0, 0.0});
    Index ids[] = {0, 1, 2};
    mesh.add_cell(ElementType::Triangle2D, ids, 3);

    auto elem = create_element(ElementType::Triangle2D);
    Assembler asm_(mesh, *elem);

    COOMatrix K;
    Vector F;
    asm_.assemble(unit_stiffness, unit_load, K, F);

    // K 是单位矩阵装配 (1 单元, 3 节点 → 对角线各 1.0)
    auto csr = coo_to_csr(K);

    Real x[] = {1.0, 2.0, 3.0};
    Real y[3];
    csr.matvec(x, y);
    EXPECT_NEAR(y[0], 1.0, 1e-15);
    EXPECT_NEAR(y[1], 2.0, 1e-15);
    EXPECT_NEAR(y[2], 3.0, 1e-15);

    // F = [1, 1, 1]
    EXPECT_NEAR(F[0], 1.0, 1e-15);
    EXPECT_NEAR(F[1], 1.0, 1e-15);
    EXPECT_NEAR(F[2], 1.0, 1e-15);
}

// ── BoundaryCondition ──
TEST(BCTest, Dirichlet_Simple) {
    // 2x2 对称矩阵 + Dirichlet BC 在节点 0: u[0]=5
    //   K = [[2, -1], [-1, 2]], F = [1, 1]
    //   施加 u[0]=5 后:
    //     行0: [1, 0], F[0]=5
    //     行1: F[1] = 1 - (-1)*5 = 6, K[1,0]=0
    COOMatrix coo;
    coo.rows = 2; coo.cols = 2;
    coo.add(0, 0,  2.0);
    coo.add(0, 1, -1.0);
    coo.add(1, 0, -1.0);
    coo.add(1, 1,  2.0);

    auto K = coo_to_csr(coo);
    Vector F = {1.0, 1.0};

    Mesh mesh;
    mesh.add_node({0.0, 0.0, 0.0});
    mesh.add_node({1.0, 0.0, 0.0});
    mesh.add_boundary("fix", {0});

    BoundaryCondition bc{BCType::Dirichlet, "fix", 5.0};
    apply_dirichlet(K, F, mesh, bc);

    // 验证行 0: 对角线=1, 其他=0, F[0]=5
    EXPECT_NEAR(F[0], 5.0, 1e-15);

    // 验证行 1: F[1] = 1 - (-1)*5 = 6
    EXPECT_NEAR(F[1], 6.0, 1e-15);
}
