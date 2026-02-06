#pragma once

#include "mesh/element.h"
#include "mesh/material.h"
#include <string>
#include <unordered_map>
#include <array>

namespace fem {

/**
 * Mesh: 单一材料域
 * 
 * 一个 Mesh 代表一块具有统一材料属性的区域。
 * 包含该区域的所有节点和单元。
 */
class Mesh {
public:
    Mesh(const std::string& name, const Material* material)
        : name_(name), material_(material), dimension_(0) {
        elements_by_dim_.fill({});
    }
    
    // ═══ 基本信息 ═══
    const std::string& name() const { return name_; }
    const Material* material() const { return material_; }
    
    int dimension() const { return dimension_; }
    
    // ═══ 节点管理 ═══
    Index add_node(const Vec3& coords) {
        Index id = nodes_.size();
        nodes_.emplace_back(coords);
        nodes_.back().set_local_id(id);
        return id;
    }
    
    const Node& node(Index id) const { return nodes_[id]; }
    Node& node(Index id) { return nodes_[id]; }
    
    std::size_t num_nodes() const { return nodes_.size(); }
    const std::vector<Node>& all_nodes() const { return nodes_; }
    
    // ═══ 单元管理 ═══
    template<typename ElemT, typename... Args>
    Index add_element(Args&&... args) {
        auto elem = std::make_unique<ElemT>(std::forward<Args>(args)...);
        
        // 确定网格维度
        int elem_dim = elem->dimension();
        if (dimension_ < elem_dim) {
            dimension_ = elem_dim;
        }
        
        // 添加到总列表
        Index id = elements_.size();
        elem->set_local_id(id);
        elements_.push_back(std::move(elem));
        
        // 按维度分类索引
        elements_by_dim_[elem_dim].push_back(id);
        
        return id;
    }
    
    const Element& element(Index id) const { return *elements_[id]; }
    Element& element(Index id) { return *elements_[id]; }
    
    std::size_t num_elements() const { return elements_.size(); }
    
    // 按维度访问
    std::size_t num_volumes() const { return elements_by_dim_[3].size(); }
    std::size_t num_faces() const { return elements_by_dim_[2].size(); }
    std::size_t num_edges() const { return elements_by_dim_[1].size(); }
    
    const std::vector<Index>& volumes() const { return elements_by_dim_[3]; }
    const std::vector<Index>& faces() const { return elements_by_dim_[2]; }
    const std::vector<Index>& edges() const { return elements_by_dim_[1]; }
    
    // ═══ 边界管理 ═══
    void add_boundary(const std::string& name, const std::vector<Index>& elem_ids) {
        boundaries_[name] = elem_ids;
    }
    
    const std::vector<Index>& boundary(const std::string& name) const {
        static std::vector<Index> empty;
        auto it = boundaries_.find(name);
        return (it != boundaries_.end()) ? it->second : empty;
    }
    
    bool has_boundary(const std::string& name) const {
        return boundaries_.find(name) != boundaries_.end();
    }
    
    std::vector<std::string> boundary_names() const {
        std::vector<std::string> names;
        for (const auto& [name, _] : boundaries_) {
            names.push_back(name);
        }
        return names;
    }
    
    // ═══ 子区域管理 ═══
    void add_region(const std::string& name, const std::vector<Index>& elem_ids) {
        regions_[name] = elem_ids;
        // 给单元打标签
        int tag = static_cast<int>(regions_.size());
        for (Index id : elem_ids) {
            elements_[id]->set_tag(tag);
        }
    }
    
    const std::vector<Index>& region(const std::string& name) const {
        static std::vector<Index> empty;
        auto it = regions_.find(name);
        return (it != regions_.end()) ? it->second : empty;
    }
    
    bool has_region(const std::string& name) const {
        return regions_.find(name) != regions_.end();
    }
    
    // ═══ 拓扑操作 ═══
    // 提取外边界
    std::vector<Index> extract_external_boundary() const;
    
    // 查找节点相邻的所有单元
    std::vector<Index> elements_containing_node(Index node_id) const;
    
    // ═══ 辅助功能 ═══
    void print_info() const;
    void validate() const;
    
private:
    std::string name_;
    const Material* material_;
    int dimension_;  // 网格的主维度 (2D or 3D)
    
    // 存储
    std::vector<Node> nodes_;
    std::vector<std::unique_ptr<Element>> elements_;
    
    // 索引
    std::array<std::vector<Index>, 4> elements_by_dim_;  // [0]=unused, [1]=edges, [2]=faces, [3]=volumes
    
    // 边界和区域
    std::unordered_map<std::string, std::vector<Index>> boundaries_;
    std::unordered_map<std::string, std::vector<Index>> regions_;
};

}  // namespace fem
