#include <gtest/gtest.h>
#include "mesh/model.h"
#include "mesh/mesh.h"
#include "mesh/material.h"
#include "mesh/element.h"

using namespace fem;

// ═══ Material Tests ═══
TEST(MaterialTest, BasicProperties) {
    Material mat(0, "Steel");
    
    mat.set_elastic(2.1e11, 0.3);
    EXPECT_DOUBLE_EQ(mat.E(), 2.1e11);
    EXPECT_DOUBLE_EQ(mat.nu(), 0.3);
    
    mat.set_thermal(50, 7850, 460);
    EXPECT_DOUBLE_EQ(mat.k(), 50);
    EXPECT_DOUBLE_EQ(mat.rho(), 7850);
    EXPECT_DOUBLE_EQ(mat.cp(), 460);
}

TEST(MaterialTest, CustomProperty) {
    Material mat(1, "Custom");
    
    mat.set_property("G", 80e9);
    EXPECT_DOUBLE_EQ(mat.property("G"), 80e9);
    
    EXPECT_DOUBLE_EQ(mat.property("nonexistent", 999), 999);
    EXPECT_FALSE(mat.has_property("nonexistent"));
    EXPECT_TRUE(mat.has_property("G"));
}

// ═══ Element Tests ═══
TEST(ElementTest, Node) {
    Node node({1.0, 2.0, 3.0});
    
    EXPECT_EQ(node.type(), ElementType::Node);
    EXPECT_EQ(node.dimension(), 0);
    EXPECT_EQ(node.num_nodes(), 1);
    EXPECT_EQ(node.order(), 1);
    
    EXPECT_DOUBLE_EQ(node.coords()[0], 1.0);
    EXPECT_DOUBLE_EQ(node.coords()[1], 2.0);
    EXPECT_DOUBLE_EQ(node.coords()[2], 3.0);
}

TEST(ElementTest, Edge2) {
    Edge2 edge(0, 1);
    
    EXPECT_EQ(edge.type(), ElementType::Edge2);
    EXPECT_EQ(edge.dimension(), 1);
    EXPECT_EQ(edge.num_nodes(), 2);
    EXPECT_EQ(edge.order(), 1);
    
    const auto& nodes = edge.nodes();
    EXPECT_EQ(nodes.size(), 2u);
    EXPECT_EQ(nodes[0], 0u);
    EXPECT_EQ(nodes[1], 1u);
}

TEST(ElementTest, Tri3) {
    Tri3 tri(0, 1, 2);
    
    EXPECT_EQ(tri.type(), ElementType::Tri3);
    EXPECT_EQ(tri.dimension(), 2);
    EXPECT_EQ(tri.num_nodes(), 3);
    EXPECT_EQ(tri.order(), 1);
    
    const auto& nodes = tri.nodes();
    EXPECT_EQ(nodes.size(), 3u);
    
    auto edges = tri.boundary_edges();
    EXPECT_EQ(edges.size(), 3u);
    EXPECT_EQ(edges[0][0], 0u);
    EXPECT_EQ(edges[0][1], 1u);
}

TEST(ElementTest, Quad4) {
    Quad4 quad(0, 1, 2, 3);
    
    EXPECT_EQ(quad.type(), ElementType::Quad4);
    EXPECT_EQ(quad.dimension(), 2);
    EXPECT_EQ(quad.num_nodes(), 4);
    EXPECT_EQ(quad.order(), 1);
    
    auto edges = quad.boundary_edges();
    EXPECT_EQ(edges.size(), 4u);
}

TEST(ElementTest, Tet4) {
    Tet4 tet(0, 1, 2, 3);
    
    EXPECT_EQ(tet.type(), ElementType::Tet4);
    EXPECT_EQ(tet.dimension(), 3);
    EXPECT_EQ(tet.num_nodes(), 4);
    EXPECT_EQ(tet.order(), 1);
    
    auto faces = tet.boundary_faces();
    EXPECT_EQ(faces.size(), 4u);
    EXPECT_EQ(faces[0].size(), 3u);
}

TEST(ElementTest, Brick8) {
    std::array<Index,8> nodes = {0,1,2,3,4,5,6,7};
    Brick8 brick(nodes);
    
    EXPECT_EQ(brick.type(), ElementType::Brick8);
    EXPECT_EQ(brick.dimension(), 3);
    EXPECT_EQ(brick.num_nodes(), 8);
    EXPECT_EQ(brick.order(), 1);
    
    auto faces = brick.boundary_faces();
    EXPECT_EQ(faces.size(), 6u);
    EXPECT_EQ(faces[0].size(), 4u);
}

// ═══ Mesh Tests ═══
TEST(MeshV2Test, BasicConstruction) {
    Material mat(0, "Steel");
    Mesh mesh("test_mesh", &mat);
    
    EXPECT_EQ(mesh.name(), "test_mesh");
    EXPECT_EQ(mesh.material(), &mat);
    EXPECT_EQ(mesh.dimension(), 0);
    EXPECT_EQ(mesh.num_nodes(), 0u);
    EXPECT_EQ(mesh.num_elements(), 0u);
}

TEST(MeshV2Test, AddNodes) {
    Material mat(0, "Steel");
    Mesh mesh("test", &mat);
    
    Index n0 = mesh.add_node({0, 0, 0});
    Index n1 = mesh.add_node({1, 0, 0});
    Index n2 = mesh.add_node({0, 1, 0});
    
    EXPECT_EQ(mesh.num_nodes(), 3u);
    EXPECT_EQ(n0, 0u);
    EXPECT_EQ(n1, 1u);
    EXPECT_EQ(n2, 2u);
    
    EXPECT_DOUBLE_EQ(mesh.node(1).coords()[0], 1.0);
}

TEST(MeshV2Test, AddElements2D) {
    Material mat(0, "Steel");
    Mesh mesh("test", &mat);
    
    Index n0 = mesh.add_node({0, 0, 0});
    Index n1 = mesh.add_node({1, 0, 0});
    Index n2 = mesh.add_node({0, 1, 0});
    
    Index f0 = mesh.add_element<Tri3>(n0, n1, n2);
    
    EXPECT_EQ(mesh.num_elements(), 1u);
    EXPECT_EQ(mesh.num_faces(), 1u);
    EXPECT_EQ(mesh.dimension(), 2);
    
    const Element& elem = mesh.element(f0);
    EXPECT_EQ(elem.type(), ElementType::Tri3);
    EXPECT_EQ(elem.local_id(), 0u);
}

TEST(MeshV2Test, AddElements3D) {
    Material mat(0, "Steel");
    Mesh mesh("test", &mat);
    
    Index n0 = mesh.add_node({0, 0, 0});
    Index n1 = mesh.add_node({1, 0, 0});
    Index n2 = mesh.add_node({0, 1, 0});
    Index n3 = mesh.add_node({0, 0, 1});
    
    Index v0 = mesh.add_element<Tet4>(n0, n1, n2, n3);
    
    EXPECT_EQ(mesh.num_elements(), 1u);
    EXPECT_EQ(mesh.num_volumes(), 1u);
    EXPECT_EQ(mesh.dimension(), 3);
}

TEST(MeshV2Test, Boundary) {
    Material mat(0, "Steel");
    Mesh mesh("test", &mat);
    
    Index n0 = mesh.add_node({0, 0, 0});
    Index n1 = mesh.add_node({1, 0, 0});
    
    Index e0 = mesh.add_element<Edge2>(n0, n1);
    
    mesh.add_boundary("left", {e0});
    
    EXPECT_TRUE(mesh.has_boundary("left"));
    EXPECT_FALSE(mesh.has_boundary("right"));
    
    const auto& boundary = mesh.boundary("left");
    EXPECT_EQ(boundary.size(), 1u);
    EXPECT_EQ(boundary[0], e0);
}

// ═══ Model Tests ═══
TEST(ModelTest, MaterialLibrary) {
    Model model("test_model");
    
    int steel_id = model.add_material("Steel");
    int aluminum_id = model.add_material("Aluminum");
    
    EXPECT_EQ(model.num_materials(), 2u);
    EXPECT_EQ(steel_id, 0);
    EXPECT_EQ(aluminum_id, 1);
    
    model.material(steel_id).set_elastic(2.1e11, 0.3);
    model.material(aluminum_id).set_elastic(7e10, 0.33);
    
    EXPECT_DOUBLE_EQ(model.material(steel_id).E(), 2.1e11);
    EXPECT_DOUBLE_EQ(model.material(aluminum_id).E(), 7e10);
    
    EXPECT_EQ(model.find_material("Steel"), steel_id);
    EXPECT_EQ(model.find_material("Aluminum"), aluminum_id);
    EXPECT_EQ(model.find_material("Nonexistent"), -1);
}

TEST(ModelTest, MultipleMeshes) {
    Model model("test_model");
    
    int steel_id = model.add_material("Steel");
    int aluminum_id = model.add_material("Aluminum");
    
    int mesh1 = model.add_mesh("steel_part", steel_id);
    int mesh2 = model.add_mesh("aluminum_part", aluminum_id);
    
    EXPECT_EQ(model.num_meshes(), 2u);
    EXPECT_EQ(mesh1, 0);
    EXPECT_EQ(mesh2, 1);
    
    EXPECT_EQ(model.mesh(mesh1).name(), "steel_part");
    EXPECT_EQ(model.mesh(mesh2).name(), "aluminum_part");
    
    EXPECT_EQ(model.mesh(mesh1).material()->name(), "Steel");
    EXPECT_EQ(model.mesh(mesh2).material()->name(), "Aluminum");
}

TEST(ModelTest, CompleteWorkflow) {
    Model model("Multi-material Structure");
    
    // 1. 定义材料
    int steel = model.add_material("Steel");
    model.material(steel).set_elastic(2.1e11, 0.3);
    model.material(steel).set_property("rho", 7850);
    
    int aluminum = model.add_material("Aluminum");
    model.material(aluminum).set_elastic(7e10, 0.33);
    model.material(aluminum).set_property("rho", 2700);
    
    // 2. 创建钢材 Mesh
    int steel_mesh_id = model.add_mesh("steel_block", steel);
    Mesh& steel_mesh = model.mesh(steel_mesh_id);
    
    Index n0 = steel_mesh.add_node({0, 0, 0});
    Index n1 = steel_mesh.add_node({1, 0, 0});
    Index n2 = steel_mesh.add_node({0, 1, 0});
    steel_mesh.add_element<Tri3>(n0, n1, n2);
    
    // 3. 创建铝材 Mesh
    int aluminum_mesh_id = model.add_mesh("aluminum_block", aluminum);
    Mesh& aluminum_mesh = model.mesh(aluminum_mesh_id);
    
    Index a0 = aluminum_mesh.add_node({1, 0, 0});
    Index a1 = aluminum_mesh.add_node({2, 0, 0});
    Index a2 = aluminum_mesh.add_node({1, 1, 0});
    aluminum_mesh.add_element<Tri3>(a0, a1, a2);
    
    // 4. 添加接触界面
    model.add_interface("steel_aluminum_contact", steel_mesh_id, aluminum_mesh_id);
    
    // 5. 验证
    EXPECT_EQ(model.total_nodes(), 6u);
    EXPECT_EQ(model.total_elements(), 2u);
    EXPECT_EQ(model.num_interfaces(), 1u);
    
    // 打印信息 (会输出到日志)
    model.print_info();
    model.validate();
}
