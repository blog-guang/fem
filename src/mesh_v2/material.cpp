#include "mesh_v2/material.h"
#include "core/logger.h"

namespace fem {
namespace v2 {

void Material::print() const {
    FEM_INFO("Material: " + name_ + " (id=" + std::to_string(id_) + ")");
    for (const auto& [key, value] : properties_) {
        FEM_INFO("  " + key + " = " + fmt_sci(value));
    }
}

}  // namespace v2
}  // namespace fem
