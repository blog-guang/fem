#include <gtest/gtest.h>
#include "mesh/mesh.h"
#include "mesh/mesh_generator.h"

using namespace fem;

// ── Mesh 基本操作 ──
TEST(MeshTest, AddNodesAndCells) {
    Mesh mesh;
    Index n0 = mesh.add_node({0.0, 0.0, 0.0});
    Index n1 = mesh.add_node({1.0, 0.0, 0.0});
    Index n2 = mesh.add_node({0.0, 1.0, 0.0});

    EXPECT_EQ(mesh.num_nodes(), 3);
    EXPECT_EQ(n0, 0); EXPECT_EQ(n1, 1); EXPECT_EQ(n2, 2);

    Index ids[] = {n0, n1, n2};
    Index c0 = mesh.add_cell(ElementType::Triangle2D, ids, 3);

    EXPECT_EQ(mesh.num_cells(), 1);
    EXPECT_EQ(c0, 0);
    EXPECT_EQ(mesh.cell(0).type, ElementType::Triangle2D);
    EXPECT_EQ(mesh.cell(0).num_nodes, 3);
}

TEST(MeshTest, CoordsAccess) {
    Mesh mesh;
    mesh.add_node({1.5, 2.5, 3.5});

    EXPECT_DOUBLE_EQ(mesh.coords(0)[0], 1.5);
    EXPECT_DOUBLE_EQ(mesh.coords(0)[1], 2.5);
    EXPECT_DOUBLE_EQ(mesh.coords(0)[2], 3.5);
}

TEST(MeshTest, Boundary) {
    Mesh mesh;
    Index n0 = mesh.add_node({0.0, 0.0, 0.0});
    Index n1 = mesh.add_node({1.0, 0.0, 0.0});
    mesh.add_boundary("bottom", {n0, n1});

    EXPECT_TRUE(mesh.has_boundary("bottom"));
    EXPECT_FALSE(mesh.has_boundary("top"));

    const auto& bnodes = mesh.boundary_nodes("bottom");
    EXPECT_EQ(bnodes.size(), 2);
    EXPECT_EQ(bnodes[0], n0);
    EXPECT_EQ(bnodes[1], n1);
}

TEST(MeshTest, BoundaryNotFound) {
    Mesh mesh;
    EXPECT_THROW(mesh.boundary_nodes("xxx"), std::invalid_argument);
}

// ── MeshGenerator ──
TEST(MeshGenTest, UnitSquareTri_1x1) {
    auto mesh = generate_unit_square_tri(1, 1);

    // 1x1: (1+1)*(1+1) = 4 节点, 1*1*2 = 2 单元
    EXPECT_EQ(mesh.num_nodes(), 4);
    EXPECT_EQ(mesh.num_cells(), 2);

    // 每个单元都是三角形
    for (std::size_t i = 0; i < mesh.num_cells(); ++i) {
        EXPECT_EQ(mesh.cell(i).type, ElementType::Triangle2D);
        EXPECT_EQ(mesh.cell(i).num_nodes, 3);
    }
}

TEST(MeshGenTest, UnitSquareTri_3x2) {
    auto mesh = generate_unit_square_tri(3, 2);

    // (3+1)*(2+1) = 12 节点, 3*2*2 = 12 单元
    EXPECT_EQ(mesh.num_nodes(), 12);
    EXPECT_EQ(mesh.num_cells(), 12);

    // 边界
    EXPECT_EQ(mesh.boundary_nodes("bottom").size(), 4);  // nx+1
    EXPECT_EQ(mesh.boundary_nodes("top").size(),    4);
    EXPECT_EQ(mesh.boundary_nodes("left").size(),   3);  // ny+1
    EXPECT_EQ(mesh.boundary_nodes("right").size(),  3);
}

TEST(MeshGenTest, UnitSquareQuad_2x2) {
    auto mesh = generate_unit_square_quad(2, 2);

    // (2+1)*(2+1) = 9 节点, 2*2 = 4 单元
    EXPECT_EQ(mesh.num_nodes(), 9);
    EXPECT_EQ(mesh.num_cells(), 4);

    for (std::size_t i = 0; i < mesh.num_cells(); ++i) {
        EXPECT_EQ(mesh.cell(i).type, ElementType::Quad2D);
        EXPECT_EQ(mesh.cell(i).num_nodes, 4);
    }
}

TEST(MeshGenTest, UnitSquareTri_CoordRange) {
    auto mesh = generate_unit_square_tri(4, 4);

    // 所有节点坐标在 [0,1] 范围内
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        EXPECT_GE(mesh.coords(i)[0], 0.0);
        EXPECT_LE(mesh.coords(i)[0], 1.0);
        EXPECT_GE(mesh.coords(i)[1], 0.0);
        EXPECT_LE(mesh.coords(i)[1], 1.0);
        EXPECT_DOUBLE_EQ(mesh.coords(i)[2], 0.0);  // 2D, z=0
    }
}
