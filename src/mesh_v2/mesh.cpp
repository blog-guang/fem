#include "mesh_v2/mesh.h"
#include "core/logger.h"
#include <set>
#include <map>
#include <algorithm>

namespace fem {
namespace v2 {

std::vector<Index> Mesh::extract_external_boundary() const {
    std::vector<Index> boundary_elems;
    
    if (dimension_ == 2) {
        // 2D: 提取边界边
        // 统计每条边出现的次数，出现1次的是边界边
        std::map<std::pair<Index,Index>, int> edge_count;
        
        for (Index fid : faces()) {
            const Face* face = dynamic_cast<const Face*>(elements_[fid].get());
            if (!face) continue;
            
            auto edges = face->boundary_edges();
            for (const auto& edge : edges) {
                Index n0 = edge[0];
                Index n1 = edge[1];
                // 规范化边的节点顺序
                if (n0 > n1) std::swap(n0, n1);
                edge_count[{n0, n1}]++;
            }
        }
        
        // 边界边只出现1次
        for (const auto& [edge, count] : edge_count) {
            if (count == 1) {
                // 创建 Edge2 单元（暂时简化，直接返回节点对）
                // 实际应用中可能需要显式创建 Edge 单元
            }
        }
        
    } else if (dimension_ == 3) {
        // 3D: 提取边界面
        // 统计每个面出现的次数，出现1次的是边界面
        // 面用节点集合表示（排序后比较）
        std::map<std::vector<Index>, int> face_count;
        
        for (Index vid : volumes()) {
            const Volume* vol = dynamic_cast<const Volume*>(elements_[vid].get());
            if (!vol) continue;
            
            auto faces = vol->boundary_faces();
            for (auto face_nodes : faces) {
                // 规范化面的节点顺序
                std::sort(face_nodes.begin(), face_nodes.end());
                face_count[face_nodes]++;
            }
        }
        
        // 边界面只出现1次
        for (const auto& [face_nodes, count] : face_count) {
            if (count == 1) {
                // 创建 Face 单元
                // 实际应用中需要根据节点数判断 Tri3/Quad4 等
            }
        }
    }
    
    return boundary_elems;
}

std::vector<Index> Mesh::elements_containing_node(Index node_id) const {
    std::vector<Index> result;
    
    for (std::size_t i = 0; i < elements_.size(); ++i) {
        const auto& nodes = elements_[i]->nodes();
        if (std::find(nodes.begin(), nodes.end(), node_id) != nodes.end()) {
            result.push_back(i);
        }
    }
    
    return result;
}

void Mesh::print_info() const {
    FEM_INFO("Mesh: " + name_);
    FEM_INFO("  Material: " + (material_ ? material_->name() : "none"));
    FEM_INFO("  Dimension: " + std::to_string(dimension_));
    FEM_INFO("  Nodes: " + std::to_string(num_nodes()));
    FEM_INFO("  Elements: " + std::to_string(num_elements()));
    if (num_volumes() > 0) FEM_INFO("    Volumes: " + std::to_string(num_volumes()));
    if (num_faces() > 0)   FEM_INFO("    Faces: " + std::to_string(num_faces()));
    if (num_edges() > 0)   FEM_INFO("    Edges: " + std::to_string(num_edges()));
    FEM_INFO("  Boundaries: " + std::to_string(boundaries_.size()));
    FEM_INFO("  Regions: " + std::to_string(regions_.size()));
}

void Mesh::validate() const {
    // 检查拓扑一致性
    bool valid = true;
    
    // 检查所有单元引用的节点是否存在
    for (std::size_t i = 0; i < elements_.size(); ++i) {
        const auto& nodes = elements_[i]->nodes();
        for (Index nid : nodes) {
            if (nid >= nodes_.size()) {
                FEM_ERROR("Element " + std::to_string(i) + " references invalid node " + std::to_string(nid));
                valid = false;
            }
        }
    }
    
    if (valid) {
        FEM_INFO("Mesh validation passed.");
    } else {
        FEM_ERROR("Mesh validation FAILED!");
    }
}

}  // namespace v2
}  // namespace fem
