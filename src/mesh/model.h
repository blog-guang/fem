#pragma once

#include "mesh/mesh.h"
#include "mesh/material.h"
#include "core/logger.h"
#include <vector>
#include <string>

namespace fem {

/**
 * Model: 顶层模型容器
 * 
 * 管理多个材料和多个 Mesh (材料域)。
 * 一个完整的 FEM 模型可能包含多个不同材料的区域。
 */
class Model {
public:
    Model(const std::string& name) : name_(name) {}
    
    const std::string& name() const { return name_; }
    
    // ═══ 材料库管理 ═══
    int add_material(const std::string& name) {
        int id = static_cast<int>(materials_.size());
        materials_.emplace_back(id, name);
        return id;
    }
    
    Material& material(int id) { return materials_[id]; }
    const Material& material(int id) const { return materials_[id]; }
    
    std::size_t num_materials() const { return materials_.size(); }
    
    int find_material(const std::string& name) const {
        for (std::size_t i = 0; i < materials_.size(); ++i) {
            if (materials_[i].name() == name) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }
    
    // ═══ Mesh 管理 ═══
    int add_mesh(const std::string& mesh_name, int material_id) {
        if (material_id < 0 || material_id >= static_cast<int>(materials_.size())) {
            FEM_ERROR("Invalid material_id: " + std::to_string(material_id));
            return -1;
        }
        
        int id = static_cast<int>(meshes_.size());
        meshes_.emplace_back(mesh_name, &materials_[material_id]);
        return id;
    }
    
    Mesh& mesh(int id) { return meshes_[id]; }
    const Mesh& mesh(int id) const { return meshes_[id]; }
    
    std::size_t num_meshes() const { return meshes_.size(); }
    
    int find_mesh(const std::string& name) const {
        for (std::size_t i = 0; i < meshes_.size(); ++i) {
            if (meshes_[i].name() == name) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }
    
    // ═══ Interface 管理 ═══
    struct Interface {
        std::string name;
        int mesh_id_1;
        int mesh_id_2;
        std::vector<std::pair<Index, Index>> node_pairs;  // 节点对应关系
        std::vector<std::pair<Index, Index>> face_pairs;  // 面对应关系
    };
    
    void add_interface(const std::string& name, int mesh1, int mesh2) {
        interfaces_.push_back({name, mesh1, mesh2, {}, {}});
    }
    
    Interface& interface(int id) { return interfaces_[id]; }
    const Interface& interface(int id) const { return interfaces_[id]; }
    
    std::size_t num_interfaces() const { return interfaces_.size(); }
    
    // ═══ 全局信息 ═══
    std::size_t total_nodes() const {
        std::size_t count = 0;
        for (const auto& mesh : meshes_) {
            count += mesh.num_nodes();
        }
        return count;
    }
    
    std::size_t total_elements() const {
        std::size_t count = 0;
        for (const auto& mesh : meshes_) {
            count += mesh.num_elements();
        }
        return count;
    }
    
    void print_info() const;
    void validate() const;
    
private:
    std::string name_;
    std::vector<Material> materials_;
    std::vector<Mesh> meshes_;
    std::vector<Interface> interfaces_;
};

}  // namespace fem
