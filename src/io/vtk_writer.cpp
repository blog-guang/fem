#include "io/vtk_writer.h"
#include "core/logger.h"
#include <stdexcept>

namespace fem {

VTKWriter::VTKWriter(const std::string& filename)
    : filename_(filename), mesh_written_(false), 
      point_data_started_(false), cell_data_started_(false),
      num_points_(0), num_cells_(0) {
    
    // 确保文件名有 .vtk 后缀
    if (filename_.size() < 4 || filename_.substr(filename_.size() - 4) != ".vtk") {
        filename_ += ".vtk";
    }
    
    file_.open(filename_);
    if (!file_.is_open()) {
        throw std::runtime_error("Failed to open VTK file: " + filename_);
    }
    
    FEM_INFO("VTK file opened: " + filename_);
}

VTKWriter::~VTKWriter() {
    close();
}

void VTKWriter::write_header(const std::string& description) {
    file_ << "# vtk DataFile Version 3.0\n";
    file_ << description << "\n";
    file_ << "ASCII\n";
    file_ << "DATASET UNSTRUCTURED_GRID\n";
}

void VTKWriter::write_points(const Mesh& mesh) {
    num_points_ = mesh.num_nodes();
    
    file_ << "POINTS " << num_points_ << " float\n";
    
    for (std::size_t i = 0; i < num_points_; ++i) {
        const auto& coords = mesh.node(i).coords();
        file_ << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
}

void VTKWriter::write_cells(const Mesh& mesh) {
    std::size_t num_cells = mesh.num_elements();
    
    // 计算 size: 每个单元的大小 = (节点数 + 1) * 单元数
    std::size_t size = 0;
    for (std::size_t i = 0; i < num_cells; ++i) {
        const auto& elem = mesh.element(i);
        size += 1 + elem.nodes().size();
    }
    
    file_ << "CELLS " << num_cells << " " << size << "\n";
    
    for (std::size_t i = 0; i < num_cells; ++i) {
        const auto& elem = mesh.element(i);
        const auto& nodes = elem.nodes();
        
        file_ << nodes.size();
        for (Index node_id : nodes) {
            file_ << " " << node_id;
        }
        file_ << "\n";
    }
}

void VTKWriter::write_cell_types(const Mesh& mesh) {
    std::size_t num_cells = mesh.num_elements();
    
    file_ << "CELL_TYPES " << num_cells << "\n";
    
    for (std::size_t i = 0; i < num_cells; ++i) {
        const auto& elem = mesh.element(i);
        int vtk_type = element_type_to_vtk(elem.type());
        file_ << vtk_type << "\n";
    }
}

int VTKWriter::element_type_to_vtk(ElementType type) const {
    // VTK 单元类型编号
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    switch (type) {
        case ElementType::Node:     return 1;   // VTK_VERTEX
        case ElementType::Edge2:    return 3;   // VTK_LINE
        case ElementType::Tri3:     return 5;   // VTK_TRIANGLE
        case ElementType::Quad4:    return 9;   // VTK_QUAD
        case ElementType::Tet4:     return 10;  // VTK_TETRA
        case ElementType::Brick8:   return 12;  // VTK_HEXAHEDRON
        default:
            FEM_WARN("Unknown element type, using VTK_VERTEX");
            return 1;
    }
}

void VTKWriter::write_mesh(const Mesh& mesh) {
    if (mesh_written_) {
        throw std::runtime_error("Mesh already written to VTK file");
    }
    
    num_cells_ = mesh.num_elements();
    
    write_header("FEM Simulation Results");
    write_points(mesh);
    write_cells(mesh);
    write_cell_types(mesh);
    
    mesh_written_ = true;
    
    FEM_INFO("Mesh written: " + std::to_string(num_points_) + " nodes, " + 
             std::to_string(num_cells_) + " elements");
}

void VTKWriter::start_point_data() {
    if (!mesh_written_) {
        throw std::runtime_error("Must write mesh before adding point data");
    }
    
    if (!point_data_started_) {
        file_ << "POINT_DATA " << num_points_ << "\n";
        point_data_started_ = true;
    }
}

void VTKWriter::add_point_scalar(const std::string& name, 
                                 const std::vector<Real>& data) {
    start_point_data();
    
    if (data.size() != num_points_) {
        throw std::runtime_error("Scalar data size mismatch: expected " + 
                               std::to_string(num_points_) + ", got " + 
                               std::to_string(data.size()));
    }
    
    file_ << "SCALARS " << name << " float 1\n";
    file_ << "LOOKUP_TABLE default\n";
    
    for (Real value : data) {
        file_ << value << "\n";
    }
    
    FEM_INFO("Added point scalar: " + name);
}

void VTKWriter::add_point_scalar(const std::string& name, 
                                 const Vector& data) {
    // 委托给 std::vector 版本
    add_point_scalar(name, data.raw());
}

void VTKWriter::add_point_vector(const std::string& name,
                                 const std::vector<Real>& data,
                                 Index dof) {
    start_point_data();
    
    if (data.size() != num_points_ * dof) {
        throw std::runtime_error("Vector data size mismatch: expected " + 
                               std::to_string(num_points_ * dof) + ", got " + 
                               std::to_string(data.size()));
    }
    
    file_ << "VECTORS " << name << " float\n";
    
    for (std::size_t i = 0; i < num_points_; ++i) {
        if (dof == 2) {
            // 2D: 补充 z=0
            file_ << data[i * 2] << " " << data[i * 2 + 1] << " 0.0\n";
        } else if (dof == 3) {
            // 3D
            file_ << data[i * 3] << " " << data[i * 3 + 1] << " " 
                  << data[i * 3 + 2] << "\n";
        } else {
            throw std::runtime_error("Unsupported DOF: " + std::to_string(dof));
        }
    }
    
    FEM_INFO("Added point vector: " + name + " (dof=" + std::to_string(dof) + ")");
}

void VTKWriter::add_point_vector(const std::string& name,
                                 const Vector& data,
                                 Index dof) {
    add_point_vector(name, data.raw(), dof);
}

void VTKWriter::start_cell_data() {
    if (!mesh_written_) {
        throw std::runtime_error("Must write mesh before adding cell data");
    }
    
    if (!cell_data_started_) {
        file_ << "CELL_DATA " << num_cells_ << "\n";
        cell_data_started_ = true;
    }
}

void VTKWriter::add_cell_scalar(const std::string& name,
                                const std::vector<Real>& data) {
    start_cell_data();
    
    if (data.size() != num_cells_) {
        throw std::runtime_error("Cell scalar data size mismatch: expected " + 
                               std::to_string(num_cells_) + ", got " + 
                               std::to_string(data.size()));
    }
    
    file_ << "SCALARS " << name << " float 1\n";
    file_ << "LOOKUP_TABLE default\n";
    
    for (Real value : data) {
        file_ << value << "\n";
    }
    
    FEM_INFO("Added cell scalar: " + name);
}

void VTKWriter::add_cell_vector(const std::string& name,
                                const std::vector<Real>& data,
                                Index dof) {
    start_cell_data();
    
    if (data.size() != num_cells_ * dof) {
        throw std::runtime_error("Cell vector data size mismatch: expected " + 
                               std::to_string(num_cells_ * dof) + ", got " + 
                               std::to_string(data.size()));
    }
    
    file_ << "VECTORS " << name << " float\n";
    
    for (std::size_t i = 0; i < num_cells_; ++i) {
        if (dof == 2) {
            // 2D: 补充 z=0
            file_ << data[i * 2] << " " << data[i * 2 + 1] << " 0.0\n";
        } else if (dof == 3) {
            // 3D
            file_ << data[i * 3] << " " << data[i * 3 + 1] << " " 
                  << data[i * 3 + 2] << "\n";
        } else {
            throw std::runtime_error("Unsupported DOF: " + std::to_string(dof));
        }
    }
    
    FEM_INFO("Added cell vector: " + name + " (dof=" + std::to_string(dof) + ")");
}

void VTKWriter::close() {
    if (file_.is_open()) {
        file_.close();
        FEM_INFO("VTK file closed: " + filename_);
    }
}

} // namespace fem
