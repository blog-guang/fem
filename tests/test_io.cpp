#include <gtest/gtest.h>
#include "io/io.h"
#include "mesh/mesh_generator.h"
#include <fstream>
#include <filesystem>

using namespace fem;
namespace fs = std::filesystem;

TEST(VTKWriterTest, WritesFile) {
    std::string dir = "/tmp/fem_test_vtk";
    auto mesh = generate_unit_square_tri(2, 2);

    // 一个标量场
    Vector temp(mesh.num_nodes(), 100.0);
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        temp[i] = static_cast<Real>(i);
    }

    VTKWriter writer(dir);
    writer.write("test.vtk", mesh, {{" temperature", &temp}});

    // 验证文件存在且非空
    std::string path = dir + "/test.vtk";
    EXPECT_TRUE(fs::exists(path));
    EXPECT_GT(fs::file_size(path), 0);

    // 验证包含关键字
    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("POINTS"), std::string::npos);
    EXPECT_NE(content.find("CELLS"), std::string::npos);
    EXPECT_NE(content.find("CELL_TYPES"), std::string::npos);
    EXPECT_NE(content.find("POINT_DATA"), std::string::npos);
    EXPECT_NE(content.find("temperature"), std::string::npos);

    // 清理
    fs::remove_all(dir);
}

TEST(VTKWriterTest, NoFields) {
    std::string dir = "/tmp/fem_test_vtk2";
    auto mesh = generate_unit_square_quad(1, 1);

    VTKWriter writer(dir);
    writer.write("no_fields.vtk", mesh);  // 不传场

    std::string path = dir + "/no_fields.vtk";
    EXPECT_TRUE(fs::exists(path));

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("POINTS"), std::string::npos);
    EXPECT_EQ(content.find("POINT_DATA"), std::string::npos);  // 无场数据

    fs::remove_all(dir);
}
