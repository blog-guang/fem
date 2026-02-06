#pragma once

#include "core/types.h"
#include <vector>
#include <string>
#include <string_view>
#include <unordered_map>

namespace fem {

// ── 单元 (轻量 POD, 无动态分配) ──
struct Cell {
    ElementType              type;
    std::array<Index, MAX_NODES> nodes{};   // 连接表
    uint8_t                  num_nodes{0};  // 实际节点数

    Index node(std::size_t i) const { return nodes[i]; }
};

// ── 网格 ──
class Mesh {
public:
    Mesh() = default;

    // ── 添加 ──
    Index add_node(const Vec3& coords);
    Index add_cell(ElementType type, const Index* node_ids, std::size_t n);

    // ── 边界标记 ──
    void add_boundary(const std::string& name, const std::vector<Index>& node_ids);
    const std::vector<Index>& boundary_nodes(std::string_view name) const;
    bool has_boundary(std::string_view name) const;

    // ── 查询 ──
    std::size_t num_nodes() const { return coords_.size(); }
    std::size_t num_cells() const { return cells_.size(); }

    const Vec3& coords(Index i) const { return coords_[i]; }
    Vec3& coords(Index i) { return coords_[i]; }
    const Cell& cell(Index i)   const { return cells_[i]; }

    const std::vector<Vec3>&  all_coords() const { return coords_; }
    const std::vector<Cell>&  all_cells()  const { return cells_; }

    void print_info() const;

private:
    std::vector<Vec3>  coords_;                          // 节点坐标
    std::vector<Cell>  cells_;                           // 单元列表
    std::unordered_map<std::string, std::vector<Index>> boundaries_;  // 边界标记
};

}  // namespace fem
