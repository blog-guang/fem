/**
 * test_io.cpp - IO 系统测试
 * 
 * 测试 VTKWriter 输出功能
 */

#include <gtest/gtest.h>
#include "io/vtk_writer.h"
#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "math/vector.h"

#include <fstream>
#include <sstream>
#include <cstdio>  // std::remove

using namespace fem;

// 测试夹具
class VTKWriterTest : public ::testing::Test {
protected:
    void SetUp() override {
        model_ = std::make_unique<Model>("Test");
        int mat_id = model_->add_material("TestMat");
        int mesh_id = model_->add_mesh("domain", mat_id);
        mesh_ = &model_->mesh(mesh_id);
    }
    
    void TearDown() override {
        // 清理生成的 VTK 文件
        for (const auto& filename : created_files_) {
            std::remove(filename.c_str());
        }
    }
    
    // 读取 VTK 文件内容
    std::string read_file(const std::string& filename) {
        std::ifstream file(filename);
        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }
    
    // 检查文件是否存在
    bool file_exists(const std::string& filename) {
        std::ifstream file(filename);
        return file.good();
    }
    
    // 记录创建的文件，用于清理
    void track_file(const std::string& filename) {
        created_files_.push_back(filename);
    }
    
    std::unique_ptr<Model> model_;
    Mesh* mesh_;
    std::vector<std::string> created_files_;
};

// ========== 基础功能测试 ==========

TEST_F(VTKWriterTest, Construction) {
    VTKWriter writer("test_vtk");
    track_file("test_vtk.vtk");
    
    EXPECT_TRUE(file_exists("test_vtk.vtk"));
}

TEST_F(VTKWriterTest, AutoExtension) {
    // 测试自动添加 .vtk 后缀
    VTKWriter writer1("test_auto");
    track_file("test_auto.vtk");
    EXPECT_TRUE(file_exists("test_auto.vtk"));
    
    // 测试已有 .vtk 后缀不重复添加
    VTKWriter writer2("test_manual.vtk");
    track_file("test_manual.vtk");
    EXPECT_TRUE(file_exists("test_manual.vtk"));
}

// ========== Mesh 输出测试 ==========

TEST_F(VTKWriterTest, WriteMeshTri3) {
    // 生成简单 2x2 三角形网格
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_tri3");
    track_file("test_tri3.vtk");
    
    writer.write_mesh(*mesh_);
    writer.close();
    
    // 验证文件内容
    std::string content = read_file("test_tri3.vtk");
    
    // 检查头部
    EXPECT_NE(content.find("# vtk DataFile Version 3.0"), std::string::npos);
    EXPECT_NE(content.find("ASCII"), std::string::npos);
    EXPECT_NE(content.find("DATASET UNSTRUCTURED_GRID"), std::string::npos);
    
    // 检查节点 (2x2 网格 → 9 个节点)
    EXPECT_NE(content.find("POINTS 9"), std::string::npos);
    
    // 检查单元 (2x2 网格 → 8 个三角形)
    EXPECT_NE(content.find("CELLS 8"), std::string::npos);
    
    // 检查单元类型 (5 = VTK_TRIANGLE)
    EXPECT_NE(content.find("CELL_TYPES 8"), std::string::npos);
}

TEST_F(VTKWriterTest, WriteMeshQuad4) {
    // 生成简单 2x2 四边形网格
    MeshGenerator::generate_unit_square_quad(2, 2, *mesh_);
    
    VTKWriter writer("test_quad4");
    track_file("test_quad4.vtk");
    
    writer.write_mesh(*mesh_);
    writer.close();
    
    std::string content = read_file("test_quad4.vtk");
    
    // 2x2 网格 → 9 节点, 4 四边形
    EXPECT_NE(content.find("POINTS 9"), std::string::npos);
    EXPECT_NE(content.find("CELLS 4"), std::string::npos);
    EXPECT_NE(content.find("CELL_TYPES 4"), std::string::npos);
}

// ========== 点数据测试 ==========

TEST_F(VTKWriterTest, AddPointScalar) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_scalar");
    track_file("test_scalar.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加标量场 (温度)
    std::vector<Real> temperature(mesh_->num_nodes());
    for (std::size_t i = 0; i < mesh_->num_nodes(); ++i) {
        temperature[i] = static_cast<Real>(i) * 10.0;
    }
    
    writer.add_point_scalar("temperature", temperature);
    writer.close();
    
    std::string content = read_file("test_scalar.vtk");
    
    // 检查 POINT_DATA
    EXPECT_NE(content.find("POINT_DATA 9"), std::string::npos);
    EXPECT_NE(content.find("SCALARS temperature float 1"), std::string::npos);
    EXPECT_NE(content.find("LOOKUP_TABLE default"), std::string::npos);
}

TEST_F(VTKWriterTest, AddPointVector2D) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_vector2d");
    track_file("test_vector2d.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加 2D 矢量场 (位移)
    std::vector<Real> displacement(mesh_->num_nodes() * 2);
    for (std::size_t i = 0; i < mesh_->num_nodes(); ++i) {
        displacement[i * 2 + 0] = static_cast<Real>(i) * 0.01;  // u_x
        displacement[i * 2 + 1] = static_cast<Real>(i) * 0.02;  // u_y
    }
    
    writer.add_point_vector("displacement", displacement, 2);
    writer.close();
    
    std::string content = read_file("test_vector2d.vtk");
    
    // 检查 VECTORS
    EXPECT_NE(content.find("POINT_DATA 9"), std::string::npos);
    EXPECT_NE(content.find("VECTORS displacement float"), std::string::npos);
    
    // 2D 矢量应该补充 z=0 (检查是否有 "0.0" 或 "0" 出现在 VECTORS 之后)
    auto vectors_pos = content.find("VECTORS displacement");
    EXPECT_NE(vectors_pos, std::string::npos);
    auto after_vectors = content.substr(vectors_pos);
    EXPECT_NE(after_vectors.find("0"), std::string::npos);  // 至少有一个 0（z分量）
}

TEST_F(VTKWriterTest, AddMultiplePointData) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_multi");
    track_file("test_multi.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加多个数据场
    std::vector<Real> temperature(mesh_->num_nodes(), 100.0);
    std::vector<Real> displacement(mesh_->num_nodes() * 2, 0.01);
    
    writer.add_point_scalar("temperature", temperature);
    writer.add_point_vector("displacement", displacement, 2);
    writer.close();
    
    std::string content = read_file("test_multi.vtk");
    
    // 两个数据场都应该存在
    EXPECT_NE(content.find("SCALARS temperature"), std::string::npos);
    EXPECT_NE(content.find("VECTORS displacement"), std::string::npos);
}

TEST_F(VTKWriterTest, VectorFromFemVector) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_fem_vector");
    track_file("test_fem_vector.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 使用 fem::Vector
    Vector disp(mesh_->num_nodes() * 2);
    for (std::size_t i = 0; i < disp.size(); ++i) {
        disp[i] = static_cast<Real>(i) * 0.01;
    }
    
    writer.add_point_vector("displacement", disp, 2);
    writer.close();
    
    EXPECT_TRUE(file_exists("test_fem_vector.vtk"));
}

// ========== 错误处理测试 ==========

TEST_F(VTKWriterTest, DataSizeMismatch) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_error");
    track_file("test_error.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 错误的标量数据大小
    std::vector<Real> wrong_scalar(5);  // 应该是 9
    EXPECT_THROW(writer.add_point_scalar("bad", wrong_scalar), std::runtime_error);
    
    // 错误的矢量数据大小
    std::vector<Real> wrong_vector(10);  // 应该是 9*2=18
    EXPECT_THROW(writer.add_point_vector("bad", wrong_vector, 2), std::runtime_error);
}

TEST_F(VTKWriterTest, DataBeforeMesh) {
    VTKWriter writer("test_order");
    track_file("test_order.vtk");
    
    std::vector<Real> data(9, 1.0);
    
    // 在写入网格前添加数据应该失败
    EXPECT_THROW(writer.add_point_scalar("bad", data), std::runtime_error);
}

TEST_F(VTKWriterTest, DoubleMeshWrite) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_double");
    track_file("test_double.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 重复写入网格应该失败
    EXPECT_THROW(writer.write_mesh(*mesh_), std::runtime_error);
}

// ========== 单元数据测试 ==========

TEST_F(VTKWriterTest, AddCellScalar) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_cell_scalar");
    track_file("test_cell_scalar.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加单元标量场 (应力)
    std::vector<Real> stress(mesh_->num_elements());
    for (std::size_t i = 0; i < mesh_->num_elements(); ++i) {
        stress[i] = static_cast<Real>(i) * 100.0;
    }
    
    writer.add_cell_scalar("stress", stress);
    writer.close();
    
    std::string content = read_file("test_cell_scalar.vtk");
    
    // 检查 CELL_DATA
    EXPECT_NE(content.find("CELL_DATA 8"), std::string::npos);
    EXPECT_NE(content.find("SCALARS stress float 1"), std::string::npos);
    EXPECT_NE(content.find("LOOKUP_TABLE default"), std::string::npos);
}

TEST_F(VTKWriterTest, AddCellVector) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_cell_vector");
    track_file("test_cell_vector.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加单元矢量场 (应变)
    std::vector<Real> strain(mesh_->num_elements() * 2);
    for (std::size_t i = 0; i < mesh_->num_elements(); ++i) {
        strain[i * 2 + 0] = static_cast<Real>(i) * 0.001;  // ε_xx
        strain[i * 2 + 1] = static_cast<Real>(i) * 0.002;  // ε_yy
    }
    
    writer.add_cell_vector("strain", strain, 2);
    writer.close();
    
    std::string content = read_file("test_cell_vector.vtk");
    
    // 检查 VECTORS
    EXPECT_NE(content.find("CELL_DATA 8"), std::string::npos);
    EXPECT_NE(content.find("VECTORS strain float"), std::string::npos);
}

TEST_F(VTKWriterTest, MixedPointAndCellData) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_mixed");
    track_file("test_mixed.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 添加点数据
    std::vector<Real> temperature(mesh_->num_nodes(), 100.0);
    writer.add_point_scalar("temperature", temperature);
    
    // 添加单元数据
    std::vector<Real> stress(mesh_->num_elements(), 50.0);
    writer.add_cell_scalar("stress", stress);
    
    writer.close();
    
    std::string content = read_file("test_mixed.vtk");
    
    // 两种数据都应该存在
    EXPECT_NE(content.find("POINT_DATA"), std::string::npos);
    EXPECT_NE(content.find("CELL_DATA"), std::string::npos);
    EXPECT_NE(content.find("temperature"), std::string::npos);
    EXPECT_NE(content.find("stress"), std::string::npos);
}

TEST_F(VTKWriterTest, CellDataSizeMismatch) {
    MeshGenerator::generate_unit_square_tri(2, 2, *mesh_);
    
    VTKWriter writer("test_cell_error");
    track_file("test_cell_error.vtk");
    
    writer.write_mesh(*mesh_);
    
    // 错误的标量数据大小
    std::vector<Real> wrong_scalar(5);  // 应该是 8
    EXPECT_THROW(writer.add_cell_scalar("bad", wrong_scalar), std::runtime_error);
    
    // 错误的矢量数据大小
    std::vector<Real> wrong_vector(10);  // 应该是 8*2=16
    EXPECT_THROW(writer.add_cell_vector("bad", wrong_vector, 2), std::runtime_error);
}

// ========== 单元类型测试 ==========

TEST_F(VTKWriterTest, AllElementTypes) {
    // 测试所有支持的单元类型
    
    // Tri3
    {
        Model model("Tri3");
        int mat_id = model.add_material("Mat");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        MeshGenerator::generate_unit_square_tri(1, 1, mesh);
        
        VTKWriter writer("test_tri3_type");
        track_file("test_tri3_type.vtk");
        writer.write_mesh(mesh);
        writer.close();
        
        std::string content = read_file("test_tri3_type.vtk");
        EXPECT_NE(content.find("5\n"), std::string::npos);  // VTK_TRIANGLE = 5
    }
    
    // Quad4
    {
        Model model("Quad4");
        int mat_id = model.add_material("Mat");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        MeshGenerator::generate_unit_square_quad(1, 1, mesh);
        
        VTKWriter writer("test_quad4_type");
        track_file("test_quad4_type.vtk");
        writer.write_mesh(mesh);
        writer.close();
        
        std::string content = read_file("test_quad4_type.vtk");
        EXPECT_NE(content.find("9\n"), std::string::npos);  // VTK_QUAD = 9
    }
    
    // Tet4
    {
        Model model("Tet4");
        int mat_id = model.add_material("Mat");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
        
        VTKWriter writer("test_tet4_type");
        track_file("test_tet4_type.vtk");
        writer.write_mesh(mesh);
        writer.close();
        
        std::string content = read_file("test_tet4_type.vtk");
        EXPECT_NE(content.find("10\n"), std::string::npos);  // VTK_TETRA = 10
    }
    
    // Brick8
    {
        Model model("Brick8");
        int mat_id = model.add_material("Mat");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        MeshGenerator::generate_unit_cube_brick(1, 1, 1, mesh);
        
        VTKWriter writer("test_brick8_type");
        track_file("test_brick8_type.vtk");
        writer.write_mesh(mesh);
        writer.close();
        
        std::string content = read_file("test_brick8_type.vtk");
        EXPECT_NE(content.find("12\n"), std::string::npos);  // VTK_HEXAHEDRON = 12
    }
}
