#include "mesh_v2/model.h"
#include "core/logger.h"

namespace fem {
namespace v2 {

void Model::print_info() const {
    FEM_INFO("Model: " + name_);
    FEM_INFO("  Materials: " + std::to_string(materials_.size()));
    FEM_INFO("  Meshes: " + std::to_string(meshes_.size()));
    FEM_INFO("  Interfaces: " + std::to_string(interfaces_.size()));
    FEM_INFO("  Total nodes: " + std::to_string(total_nodes()));
    FEM_INFO("  Total elements: " + std::to_string(total_elements()));
    
    FEM_INFO("");
    FEM_INFO("=== Materials ===");
    for (const auto& mat : materials_) {
        mat.print();
    }
    
    FEM_INFO("");
    FEM_INFO("=== Meshes ===");
    for (std::size_t i = 0; i < meshes_.size(); ++i) {
        FEM_INFO("--- Mesh " + std::to_string(i) + " ---");
        meshes_[i].print_info();
    }
}

void Model::validate() const {
    FEM_INFO("Validating model: " + name_);
    
    // 验证每个 Mesh
    for (std::size_t i = 0; i < meshes_.size(); ++i) {
        FEM_INFO("Validating Mesh " + std::to_string(i) + ": " + meshes_[i].name());
        meshes_[i].validate();
    }
    
    // 验证 Interface 引用的 Mesh 是否存在
    for (std::size_t i = 0; i < interfaces_.size(); ++i) {
        const auto& iface = interfaces_[i];
        if (iface.mesh_id_1 < 0 || iface.mesh_id_1 >= static_cast<int>(meshes_.size())) {
            FEM_ERROR("Interface " + std::to_string(i) + " references invalid mesh_id_1: " + std::to_string(iface.mesh_id_1));
        }
        if (iface.mesh_id_2 < 0 || iface.mesh_id_2 >= static_cast<int>(meshes_.size())) {
            FEM_ERROR("Interface " + std::to_string(i) + " references invalid mesh_id_2: " + std::to_string(iface.mesh_id_2));
        }
    }
    
    FEM_INFO("Model validation complete.");
}

}  // namespace v2
}  // namespace fem
