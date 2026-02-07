/**
 * vtk_writer.h - VTK 文件输出
 * 
 * 支持输出有限元结果到 VTK 格式，用于 ParaView 等可视化工具
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "math/vector.h"

#include <string>
#include <vector>
#include <fstream>

namespace fem {

/**
 * VTK 输出器
 * 
 * 使用方法:
 * ```cpp
 * VTKWriter writer("output.vtk");
 * writer.write_mesh(mesh);
 * writer.add_point_scalar("temperature", temperature_data);
 * writer.add_point_vector("displacement", displacement_data, 2);
 * writer.close();
 * ```
 */
class VTKWriter {
public:
    /**
     * 构造函数
     * @param filename 输出文件名 (自动添加 .vtk 后缀)
     */
    explicit VTKWriter(const std::string& filename);

    /**
     * 析构函数 (自动关闭文件)
     */
    ~VTKWriter();

    /**
     * 写入网格结构
     * @param mesh 网格对象
     */
    void write_mesh(const Mesh& mesh);

    /**
     * 添加节点标量场数据
     * @param name 数据名称
     * @param data 数据向量 (大小应等于节点数)
     */
    void add_point_scalar(const std::string& name, const std::vector<Real>& data);

    /**
     * 添加节点矢量场数据
     * @param name 数据名称
     * @param data 数据向量 (大小应为节点数 * dof)
     * @param dof 每个节点的自由度数 (2 或 3)
     */
    void add_point_vector(const std::string& name, 
                         const std::vector<Real>& data,
                         Index dof);

    /**
     * 添加节点矢量场数据 (fem::Vector 版本)
     */
    void add_point_vector(const std::string& name,
                         const Vector& data,
                         Index dof);

    /**
     * 添加单元标量场数据
     * @param name 数据名称
     * @param data 数据向量 (大小应等于单元数)
     */
    void add_cell_scalar(const std::string& name, const std::vector<Real>& data);

    /**
     * 添加单元矢量场数据
     * @param name 数据名称
     * @param data 数据向量 (大小应为单元数 * dof)
     * @param dof 每个单元的维度 (2 或 3)
     */
    void add_cell_vector(const std::string& name,
                        const std::vector<Real>& data,
                        Index dof);

    /**
     * 关闭文件
     */
    void close();

private:
    std::ofstream file_;
    std::string filename_;
    bool mesh_written_;
    bool point_data_started_;
    bool cell_data_started_;
    std::size_t num_points_;
    std::size_t num_cells_;

    /**
     * 写入 VTK 文件头
     */
    void write_header(const std::string& description);

    /**
     * 写入节点坐标
     */
    void write_points(const Mesh& mesh);

    /**
     * 写入单元连接性
     */
    void write_cells(const Mesh& mesh);

    /**
     * 写入单元类型
     */
    void write_cell_types(const Mesh& mesh);

    /**
     * 开始 POINT_DATA 部分
     */
    void start_point_data();

    /**
     * 开始 CELL_DATA 部分
     */
    void start_cell_data();

    /**
     * 将 ElementType 转换为 VTK 单元类型编号
     */
    int element_type_to_vtk(ElementType type) const;
};

} // namespace fem
