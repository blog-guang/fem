#include <gtest/gtest.h>
#include "mesh/mesh_generator.h"
#include "mesh/model.h"

using namespace fem;

// ═══ 2D Triangle Mesh Tests ═══
TEST(MeshGeneratorTest, UnitSquareTri_1x1) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(1, 1, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 4u);     // 2x2 nodes
    EXPECT_EQ(mesh.num_elements(), 2u);  // 2 triangles
    EXPECT_EQ(mesh.num_faces(), 2u);
    EXPECT_EQ(mesh.dimension(), 2);
}

TEST(MeshGeneratorTest, UnitSquareTri_3x2) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(3, 2, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 12u);    // 4x3 nodes
    EXPECT_EQ(mesh.num_elements(), 12u); // 6 cells * 2 triangles
    EXPECT_EQ(mesh.num_faces(), 12u);
}

TEST(MeshGeneratorTest, UnitSquareTri_CoordCheck) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    // 检查角点坐标
    const auto& n0 = mesh.node(0).coords();
    EXPECT_DOUBLE_EQ(n0[0], 0.0);
    EXPECT_DOUBLE_EQ(n0[1], 0.0);
    
    const auto& n2 = mesh.node(2).coords();
    EXPECT_DOUBLE_EQ(n2[0], 1.0);
    EXPECT_DOUBLE_EQ(n2[1], 0.0);
    
    const auto& n8 = mesh.node(8).coords();
    EXPECT_DOUBLE_EQ(n8[0], 1.0);
    EXPECT_DOUBLE_EQ(n8[1], 1.0);
}

// ═══ 2D Quad Mesh Tests ═══
TEST(MeshGeneratorTest, UnitSquareQuad_2x2) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_quad(2, 2, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 9u);    // 3x3 nodes
    EXPECT_EQ(mesh.num_elements(), 4u); // 4 quads
    EXPECT_EQ(mesh.num_faces(), 4u);
    EXPECT_EQ(mesh.dimension(), 2);
}

// ═══ 3D Tet Mesh Tests ═══
TEST(MeshGeneratorTest, UnitCubeTet_1x1x1) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 8u);    // 2x2x2 nodes
    EXPECT_EQ(mesh.num_elements(), 5u); // 5 tetrahedra per cube
    EXPECT_EQ(mesh.num_volumes(), 5u);
    EXPECT_EQ(mesh.dimension(), 3);
}

TEST(MeshGeneratorTest, UnitCubeTet_2x2x2) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(2, 2, 2, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 27u);   // 3x3x3 nodes
    EXPECT_EQ(mesh.num_elements(), 40u); // 8 cubes * 5 tets
    EXPECT_EQ(mesh.num_volumes(), 40u);
}

// ═══ 3D Brick Mesh Tests ═══
TEST(MeshGeneratorTest, UnitCubeBrick_1x1x1) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_brick(1, 1, 1, mesh);
    
    EXPECT_EQ(mesh.num_nodes(), 8u);    // 2x2x2 nodes
    EXPECT_EQ(mesh.num_elements(), 1u); // 1 brick
    EXPECT_EQ(mesh.num_volumes(), 1u);
    EXPECT_EQ(mesh.dimension(), 3);
}

TEST(MeshGeneratorTest, UnitCubeBrick_CoordCheck) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_brick(1, 1, 1, mesh);
    
    // 检查8个角点
    const auto& n0 = mesh.node(0).coords();
    EXPECT_DOUBLE_EQ(n0[0], 0.0);
    EXPECT_DOUBLE_EQ(n0[1], 0.0);
    EXPECT_DOUBLE_EQ(n0[2], 0.0);
    
    const auto& n7 = mesh.node(7).coords();
    EXPECT_DOUBLE_EQ(n7[0], 1.0);
    EXPECT_DOUBLE_EQ(n7[1], 1.0);
    EXPECT_DOUBLE_EQ(n7[2], 1.0);
}

// ═══ Boundary Identification Tests ═══
TEST(MeshGeneratorTest, IdentifyBoundaries2D) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    EXPECT_TRUE(mesh.has_boundary("left"));
    EXPECT_TRUE(mesh.has_boundary("right"));
    EXPECT_TRUE(mesh.has_boundary("bottom"));
    EXPECT_TRUE(mesh.has_boundary("top"));
    
    // 检查边界节点数
    EXPECT_EQ(mesh.boundary("left").size(), 3u);    // 3 nodes on left
    EXPECT_EQ(mesh.boundary("right").size(), 3u);   // 3 nodes on right
    EXPECT_EQ(mesh.boundary("bottom").size(), 3u);  // 3 nodes on bottom
    EXPECT_EQ(mesh.boundary("top").size(), 3u);     // 3 nodes on top
}

TEST(MeshGeneratorTest, IdentifyBoundaries3D) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_brick(1, 1, 1, mesh);
    MeshGenerator::identify_boundaries_3d(mesh);
    
    EXPECT_TRUE(mesh.has_boundary("left"));
    EXPECT_TRUE(mesh.has_boundary("right"));
    EXPECT_TRUE(mesh.has_boundary("bottom"));
    EXPECT_TRUE(mesh.has_boundary("top"));
    EXPECT_TRUE(mesh.has_boundary("front"));
    EXPECT_TRUE(mesh.has_boundary("back"));
    
    // 每个面4个节点
    EXPECT_EQ(mesh.boundary("left").size(), 4u);
    EXPECT_EQ(mesh.boundary("right").size(), 4u);
}

// ═══ Element Type Tests ═══
TEST(MeshGeneratorTest, ElementTypes) {
    Model model("test");
    int mat_id = model.add_material("TestMat");
    
    // Triangle mesh
    int tri_mesh_id = model.add_mesh("tri_mesh", mat_id);
    Mesh& tri_mesh = model.mesh(tri_mesh_id);
    MeshGenerator::generate_unit_square_tri(1, 1, tri_mesh);
    EXPECT_EQ(tri_mesh.element(0).type(), ElementType::Tri3);
    
    // Quad mesh
    int quad_mesh_id = model.add_mesh("quad_mesh", mat_id);
    Mesh& quad_mesh = model.mesh(quad_mesh_id);
    MeshGenerator::generate_unit_square_quad(1, 1, quad_mesh);
    EXPECT_EQ(quad_mesh.element(0).type(), ElementType::Quad4);
    
    // Tet mesh
    int tet_mesh_id = model.add_mesh("tet_mesh", mat_id);
    Mesh& tet_mesh = model.mesh(tet_mesh_id);
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, tet_mesh);
    EXPECT_EQ(tet_mesh.element(0).type(), ElementType::Tet4);
    
    // Brick mesh
    int brick_mesh_id = model.add_mesh("brick_mesh", mat_id);
    Mesh& brick_mesh = model.mesh(brick_mesh_id);
    MeshGenerator::generate_unit_cube_brick(1, 1, 1, brick_mesh);
    EXPECT_EQ(brick_mesh.element(0).type(), ElementType::Brick8);
}
