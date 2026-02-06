#include "mesh/mesh_generator.h"
#include "core/logger.h"
#include <cmath>

namespace fem {

void MeshGenerator::generate_unit_square_tri(int nx, int ny, Mesh& mesh) {
    if (nx <= 0 || ny <= 0) {
        FEM_ERROR("Invalid mesh size: nx=" + std::to_string(nx) + ", ny=" + std::to_string(ny));
        return;
    }
    
    int n_nodes_x = nx + 1;
    int n_nodes_y = ny + 1;
    Real dx = 1.0 / nx;
    Real dy = 1.0 / ny;
    
    // 添加节点
    for (int j = 0; j < n_nodes_y; ++j) {
        for (int i = 0; i < n_nodes_x; ++i) {
            Real x = i * dx;
            Real y = j * dy;
            mesh.add_node({x, y, 0.0});
        }
    }
    
    // 添加三角形单元 (每个矩形分为2个三角形)
    auto node_id = [n_nodes_x](int i, int j) -> Index {
        return j * n_nodes_x + i;
    };
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Index n00 = node_id(i, j);
            Index n10 = node_id(i+1, j);
            Index n11 = node_id(i+1, j+1);
            Index n01 = node_id(i, j+1);
            
            // 下三角形
            mesh.add_element<Tri3>(n00, n10, n11);
            // 上三角形
            mesh.add_element<Tri3>(n00, n11, n01);
        }
    }
    
    FEM_INFO("Generated unit_square_tri: " + std::to_string(nx) + "x" + std::to_string(ny));
    FEM_INFO("  Nodes: " + std::to_string(mesh.num_nodes()));
    FEM_INFO("  Elements: " + std::to_string(mesh.num_elements()));
}

void MeshGenerator::generate_unit_square_quad(int nx, int ny, Mesh& mesh) {
    if (nx <= 0 || ny <= 0) {
        FEM_ERROR("Invalid mesh size: nx=" + std::to_string(nx) + ", ny=" + std::to_string(ny));
        return;
    }
    
    int n_nodes_x = nx + 1;
    int n_nodes_y = ny + 1;
    Real dx = 1.0 / nx;
    Real dy = 1.0 / ny;
    
    // 添加节点
    for (int j = 0; j < n_nodes_y; ++j) {
        for (int i = 0; i < n_nodes_x; ++i) {
            Real x = i * dx;
            Real y = j * dy;
            mesh.add_node({x, y, 0.0});
        }
    }
    
    // 添加四边形单元
    auto node_id = [n_nodes_x](int i, int j) -> Index {
        return j * n_nodes_x + i;
    };
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Index n00 = node_id(i, j);
            Index n10 = node_id(i+1, j);
            Index n11 = node_id(i+1, j+1);
            Index n01 = node_id(i, j+1);
            
            mesh.add_element<Quad4>(n00, n10, n11, n01);
        }
    }
    
    FEM_INFO("Generated unit_square_quad: " + std::to_string(nx) + "x" + std::to_string(ny));
    FEM_INFO("  Nodes: " + std::to_string(mesh.num_nodes()));
    FEM_INFO("  Elements: " + std::to_string(mesh.num_elements()));
}

void MeshGenerator::generate_unit_cube_tet(int nx, int ny, int nz, Mesh& mesh) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        FEM_ERROR("Invalid mesh size: nx=" + std::to_string(nx) + 
                 ", ny=" + std::to_string(ny) + ", nz=" + std::to_string(nz));
        return;
    }
    
    int n_nodes_x = nx + 1;
    int n_nodes_y = ny + 1;
    int n_nodes_z = nz + 1;
    Real dx = 1.0 / nx;
    Real dy = 1.0 / ny;
    Real dz = 1.0 / nz;
    
    // 添加节点
    for (int k = 0; k < n_nodes_z; ++k) {
        for (int j = 0; j < n_nodes_y; ++j) {
            for (int i = 0; i < n_nodes_x; ++i) {
                Real x = i * dx;
                Real y = j * dy;
                Real z = k * dz;
                mesh.add_node({x, y, z});
            }
        }
    }
    
    // 添加四面体单元 (每个六面体分为5个四面体)
    auto node_id = [n_nodes_x, n_nodes_y](int i, int j, int k) -> Index {
        return k * (n_nodes_x * n_nodes_y) + j * n_nodes_x + i;
    };
    
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // 六面体的8个顶点
                Index n000 = node_id(i,   j,   k);
                Index n100 = node_id(i+1, j,   k);
                Index n010 = node_id(i,   j+1, k);
                Index n110 = node_id(i+1, j+1, k);
                Index n001 = node_id(i,   j,   k+1);
                Index n101 = node_id(i+1, j,   k+1);
                Index n011 = node_id(i,   j+1, k+1);
                Index n111 = node_id(i+1, j+1, k+1);
                
                // 5个四面体划分 (保持一致的方向)
                mesh.add_element<Tet4>(n000, n100, n110, n111);
                mesh.add_element<Tet4>(n000, n100, n111, n101);
                mesh.add_element<Tet4>(n000, n111, n011, n001);
                mesh.add_element<Tet4>(n000, n110, n011, n111);
                mesh.add_element<Tet4>(n000, n010, n110, n011);
            }
        }
    }
    
    FEM_INFO("Generated unit_cube_tet: " + std::to_string(nx) + "x" + 
             std::to_string(ny) + "x" + std::to_string(nz));
    FEM_INFO("  Nodes: " + std::to_string(mesh.num_nodes()));
    FEM_INFO("  Elements: " + std::to_string(mesh.num_elements()));
}

void MeshGenerator::generate_unit_cube_brick(int nx, int ny, int nz, Mesh& mesh) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        FEM_ERROR("Invalid mesh size: nx=" + std::to_string(nx) + 
                 ", ny=" + std::to_string(ny) + ", nz=" + std::to_string(nz));
        return;
    }
    
    int n_nodes_x = nx + 1;
    int n_nodes_y = ny + 1;
    int n_nodes_z = nz + 1;
    Real dx = 1.0 / nx;
    Real dy = 1.0 / ny;
    Real dz = 1.0 / nz;
    
    // 添加节点
    for (int k = 0; k < n_nodes_z; ++k) {
        for (int j = 0; j < n_nodes_y; ++j) {
            for (int i = 0; i < n_nodes_x; ++i) {
                Real x = i * dx;
                Real y = j * dy;
                Real z = k * dz;
                mesh.add_node({x, y, z});
            }
        }
    }
    
    // 添加六面体单元
    auto node_id = [n_nodes_x, n_nodes_y](int i, int j, int k) -> Index {
        return k * (n_nodes_x * n_nodes_y) + j * n_nodes_x + i;
    };
    
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                Index n000 = node_id(i,   j,   k);
                Index n100 = node_id(i+1, j,   k);
                Index n110 = node_id(i+1, j+1, k);
                Index n010 = node_id(i,   j+1, k);
                Index n001 = node_id(i,   j,   k+1);
                Index n101 = node_id(i+1, j,   k+1);
                Index n111 = node_id(i+1, j+1, k+1);
                Index n011 = node_id(i,   j+1, k+1);
                
                std::array<Index,8> nodes = {n000, n100, n110, n010, 
                                             n001, n101, n111, n011};
                mesh.add_element<Brick8>(nodes);
            }
        }
    }
    
    FEM_INFO("Generated unit_cube_brick: " + std::to_string(nx) + "x" + 
             std::to_string(ny) + "x" + std::to_string(nz));
    FEM_INFO("  Nodes: " + std::to_string(mesh.num_nodes()));
    FEM_INFO("  Elements: " + std::to_string(mesh.num_elements()));
}

void MeshGenerator::identify_boundaries_2d(Mesh& mesh) {
    if (mesh.num_nodes() == 0) return;
    
    Real tol = 1e-10;
    std::vector<Index> left, right, bottom, top;
    
    // 遍历所有节点，根据坐标分类
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        const auto& coords = mesh.node(i).coords();
        Real x = coords[0];
        Real y = coords[1];
        
        if (std::abs(x - 0.0) < tol) left.push_back(i);
        if (std::abs(x - 1.0) < tol) right.push_back(i);
        if (std::abs(y - 0.0) < tol) bottom.push_back(i);
        if (std::abs(y - 1.0) < tol) top.push_back(i);
    }
    
    // 为边界节点创建 Edge 单元
    auto create_boundary_edges = [&mesh](const std::vector<Index>& nodes) -> std::vector<Index> {
        std::vector<Index> edge_ids;
        // 简化: 暂时不创建 Edge 单元，只记录节点
        // 完整实现需要根据拓扑连接创建边
        return edge_ids;
    };
    
    // 暂时用节点列表作为边界（后续可改进为真正的 Edge 单元）
    if (!left.empty())   mesh.add_boundary("left", left);
    if (!right.empty())  mesh.add_boundary("right", right);
    if (!bottom.empty()) mesh.add_boundary("bottom", bottom);
    if (!top.empty())    mesh.add_boundary("top", top);
    
    FEM_INFO("Identified 2D boundaries:");
    if (!left.empty())   FEM_INFO("  left: " + std::to_string(left.size()) + " nodes");
    if (!right.empty())  FEM_INFO("  right: " + std::to_string(right.size()) + " nodes");
    if (!bottom.empty()) FEM_INFO("  bottom: " + std::to_string(bottom.size()) + " nodes");
    if (!top.empty())    FEM_INFO("  top: " + std::to_string(top.size()) + " nodes");
}

void MeshGenerator::identify_boundaries_3d(Mesh& mesh) {
    if (mesh.num_nodes() == 0) return;
    
    Real tol = 1e-10;
    std::vector<Index> left, right, bottom, top, front, back;
    
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        const auto& coords = mesh.node(i).coords();
        Real x = coords[0];
        Real y = coords[1];
        Real z = coords[2];
        
        if (std::abs(x - 0.0) < tol) left.push_back(i);
        if (std::abs(x - 1.0) < tol) right.push_back(i);
        if (std::abs(y - 0.0) < tol) bottom.push_back(i);
        if (std::abs(y - 1.0) < tol) top.push_back(i);
        if (std::abs(z - 0.0) < tol) front.push_back(i);
        if (std::abs(z - 1.0) < tol) back.push_back(i);
    }
    
    if (!left.empty())   mesh.add_boundary("left", left);
    if (!right.empty())  mesh.add_boundary("right", right);
    if (!bottom.empty()) mesh.add_boundary("bottom", bottom);
    if (!top.empty())    mesh.add_boundary("top", top);
    if (!front.empty())  mesh.add_boundary("front", front);
    if (!back.empty())   mesh.add_boundary("back", back);
    
    FEM_INFO("Identified 3D boundaries:");
    if (!left.empty())   FEM_INFO("  left: " + std::to_string(left.size()) + " nodes");
    if (!right.empty())  FEM_INFO("  right: " + std::to_string(right.size()) + " nodes");
    if (!bottom.empty()) FEM_INFO("  bottom: " + std::to_string(bottom.size()) + " nodes");
    if (!top.empty())    FEM_INFO("  top: " + std::to_string(top.size()) + " nodes");
    if (!front.empty())  FEM_INFO("  front: " + std::to_string(front.size()) + " nodes");
    if (!back.empty())   FEM_INFO("  back: " + std::to_string(back.size()) + " nodes");
}

}  // namespace fem
