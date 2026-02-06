#include "mesh/mesh.h"
#include "core/logger.h"
#include <stdexcept>

namespace fem {

Index Mesh::add_node(const Vec3& coords) {
    Index id = coords_.size();
    coords_.push_back(coords);
    return id;
}

Index Mesh::add_cell(ElementType type, const Index* node_ids, std::size_t n) {
    if (n > MAX_NODES) {
        throw std::invalid_argument("add_cell: n > MAX_NODES");
    }
    Cell c;
    c.type      = type;
    c.num_nodes = static_cast<uint8_t>(n);
    for (std::size_t i = 0; i < n; ++i) {
        c.nodes[i] = node_ids[i];
    }
    Index id = cells_.size();
    cells_.push_back(c);
    return id;
}

void Mesh::add_boundary(const std::string& name, const std::vector<Index>& node_ids) {
    boundaries_[name] = node_ids;
}

const std::vector<Index>& Mesh::boundary_nodes(std::string_view name) const {
    auto it = boundaries_.find(std::string(name));
    if (it == boundaries_.end()) {
        throw std::invalid_argument("boundary '" + std::string(name) + "' not found");
    }
    return it->second;
}

bool Mesh::has_boundary(std::string_view name) const {
    return boundaries_.count(std::string(name)) > 0;
}

void Mesh::print_info() const {
    FEM_INFO("Mesh: " + std::to_string(num_nodes()) + " nodes, " +
             std::to_string(num_cells()) + " cells");
}

}  // namespace fem
