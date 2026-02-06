#include "io/io.h"
#include "core/logger.h"
#include <fstream>
#include <stdexcept>
#include <filesystem>

namespace fem {

namespace fs = std::filesystem;

// VTK 单元类型 ID
static int vtk_cell_type(ElementType type) {
    switch (type) {
        case ElementType::Triangle2D:    return 5;   // VTK_TRIANGLE
        case ElementType::Quad2D:        return 9;   // VTK_QUAD
        case ElementType::Tetrahedron3D: return 10;  // VTK_TETRA
        case ElementType::Hexahedron3D:  return 12;  // VTK_HEXAHEDRON
    }
    return -1;
}

VTKWriter::VTKWriter(const std::string& output_dir)
    : output_dir_(output_dir)
{
    fs::create_directories(output_dir_);
}

void VTKWriter::write(const std::string& filename,
                       const Mesh&        mesh,
                       const std::vector<ScalarField>& fields) const
{
    std::string path = (fs::path(output_dir_) / filename).string();
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("VTKWriter: cannot open " + path);
    }

    std::size_t np = mesh.num_nodes();
    std::size_t nc = mesh.num_cells();

    // ── Header ──
    out << "# vtk DataFile Version 3.0\n";
    out << "fem_multiphysics\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // ── 节点 ──
    out << "POINTS " << np << " double\n";
    for (std::size_t i = 0; i < np; ++i) {
        const auto& c = mesh.coords(i);
        out << c[0] << " " << c[1] << " " << c[2] << "\n";
    }

    // ── 单元连接 ──
    std::size_t total_conn = 0;
    for (std::size_t i = 0; i < nc; ++i) {
        total_conn += 1 + mesh.cell(i).num_nodes;
    }
    out << "\nCELLS " << nc << " " << total_conn << "\n";
    for (std::size_t i = 0; i < nc; ++i) {
        const auto& c = mesh.cell(i);
        out << static_cast<int>(c.num_nodes);
        for (uint8_t j = 0; j < c.num_nodes; ++j) {
            out << " " << c.node(j);
        }
        out << "\n";
    }

    // ── 单元类型 ──
    out << "\nCELL_TYPES " << nc << "\n";
    for (std::size_t i = 0; i < nc; ++i) {
        out << vtk_cell_type(mesh.cell(i).type) << "\n";
    }

    // ── 节点数据 (标量场) ──
    if (!fields.empty()) {
        out << "\nPOINT_DATA " << np << "\n";
        for (const auto& f : fields) {
            out << "SCALARS " << f.name << " double 1\n";
            out << "LOOKUP_TABLE default\n";
            for (std::size_t i = 0; i < np; ++i) {
                out << (*f.values)[i] << "\n";
            }
        }
    }

    FEM_INFO("VTK written: " + path);
}

}  // namespace fem
