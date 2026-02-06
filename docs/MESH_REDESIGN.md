# Mesh 数据结构重构设计

## 概述

重新设计 FEM 网格数据结构，采用层次化的 Element 体系，支持：
- 拓扑关系管理 (Volume → Face → Edge → Node)
- 材料属性管理
- 多区域/多材料支持
- 高阶单元扩展性

---

## 1. Element 类层次

### 1.1 Element (抽象基类)

```cpp
class Element {
public:
    virtual ~Element() = default;
    
    // 基本属性
    virtual ElementType type() const = 0;
    virtual int dimension() const = 0;      // 0=Node, 1=Edge, 2=Face, 3=Volume
    virtual int num_nodes() const = 0;
    
    // 节点访问
    virtual const std::vector<Index>& node_ids() const = 0;
    
    // 拓扑关系
    virtual int num_edges() const { return 0; }
    virtual int num_faces() const { return 0; }
    
    // 标签/属性
    void set_tag(int tag) { tag_ = tag; }
    int tag() const { return tag_; }
    
protected:
    int tag_{0};  // 用于区域/材料标记
};
```

### 1.2 Node (0D)

```cpp
class Node : public Element {
public:
    Node(Index id, const Vec3& coords);
    
    ElementType type() const override { return ElementType::Node; }
    int dimension() const override { return 0; }
    int num_nodes() const override { return 1; }
    
    const Vec3& coords() const { return coords_; }
    Vec3& coords() { return coords_; }
    
    Index id() const { return id_; }
    
private:
    Index id_;
    Vec3 coords_;
};
```

### 1.3 Edge (1D)

```cpp
class Edge : public Element {
public:
    virtual int order() const = 0;  // 1=线性, 2=二次
};

class Edge2 : public Edge {
public:
    Edge2(Index n0, Index n1);
    
    ElementType type() const override { return ElementType::Edge2; }
    int dimension() const override { return 1; }
    int num_nodes() const override { return 2; }
    int order() const override { return 1; }
    
    const std::vector<Index>& node_ids() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;  // [n0, n1]
};

class Edge3 : public Edge {
public:
    Edge3(Index n0, Index n1, Index n_mid);
    
    ElementType type() const override { return ElementType::Edge3; }
    int num_nodes() const override { return 3; }
    int order() const override { return 2; }
    
private:
    std::vector<Index> nodes_;  // [n0, n1, n_mid]
};
```

### 1.4 Face (2D)

```cpp
class Face : public Element {
public:
    virtual int order() const = 0;
    int dimension() const override { return 2; }
    
    // 面的边界边
    virtual std::vector<std::array<Index,2>> boundary_edges() const = 0;
};

// 3节点三角形
class Tri3 : public Face {
public:
    Tri3(Index n0, Index n1, Index n2);
    
    ElementType type() const override { return ElementType::Tri3; }
    int num_nodes() const override { return 3; }
    int num_edges() const override { return 3; }
    int order() const override { return 1; }
    
    const std::vector<Index>& node_ids() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;  // [n0, n1, n2]
};

// 6节点三角形 (二次)
class Tri6 : public Face {
public:
    Tri6(const std::array<Index,6>& nodes);
    
    ElementType type() const override { return ElementType::Tri6; }
    int num_nodes() const override { return 6; }
    int num_edges() const override { return 3; }
    int order() const override { return 2; }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2, n01,n12,n20]
};

// 4节点四边形
class Quad4 : public Face {
public:
    Quad4(Index n0, Index n1, Index n2, Index n3);
    
    ElementType type() const override { return ElementType::Quad4; }
    int num_nodes() const override { return 4; }
    int num_edges() const override { return 4; }
    int order() const override { return 1; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[3]},
                {nodes_[3], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 8节点四边形 (二次)
class Quad8 : public Face {
public:
    Quad8(const std::array<Index,8>& nodes);
    
    ElementType type() const override { return ElementType::Quad8; }
    int num_nodes() const override { return 8; }
    int order() const override { return 2; }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2,n3, n01,n12,n23,n30]
};
```

### 1.5 Volume (3D)

```cpp
class Volume : public Element {
public:
    virtual int order() const = 0;
    int dimension() const override { return 3; }
    
    // 体单元的边界面
    virtual std::vector<std::vector<Index>> boundary_faces() const = 0;
    
    // 材料标记
    void set_material_id(int mat_id) { material_id_ = mat_id; }
    int material_id() const { return material_id_; }
    
protected:
    int material_id_{0};
};

// 4节点四面体
class Tet4 : public Volume {
public:
    Tet4(Index n0, Index n1, Index n2, Index n3);
    
    ElementType type() const override { return ElementType::Tet4; }
    int num_nodes() const override { return 4; }
    int num_faces() const override { return 4; }
    int order() const override { return 1; }
    
    const std::vector<Index>& node_ids() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[2], nodes_[1]},  // 底面
                {nodes_[0], nodes_[1], nodes_[3]},  // 侧面1
                {nodes_[1], nodes_[2], nodes_[3]},  // 侧面2
                {nodes_[2], nodes_[0], nodes_[3]}}; // 侧面3
    }
    
private:
    std::vector<Index> nodes_;
};

// 10节点四面体 (二次)
class Tet10 : public Volume {
public:
    Tet10(const std::array<Index,10>& nodes);
    
    ElementType type() const override { return ElementType::Tet10; }
    int num_nodes() const override { return 10; }
    int num_faces() const override { return 4; }
    int order() const override { return 2; }
    
private:
    std::vector<Index> nodes_;  // 4 corner + 6 edge
};

// 8节点六面体 (Brick)
class Brick8 : public Volume {
public:
    Brick8(const std::array<Index,8>& nodes);
    
    ElementType type() const override { return ElementType::Brick8; }
    int num_nodes() const override { return 8; }
    int num_faces() const override { return 6; }
    int order() const override { return 1; }
    
    const std::vector<Index>& node_ids() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[3], nodes_[2], nodes_[1]},  // 底面
                {nodes_[4], nodes_[5], nodes_[6], nodes_[7]},  // 顶面
                {nodes_[0], nodes_[1], nodes_[5], nodes_[4]},  // 前面
                {nodes_[2], nodes_[3], nodes_[7], nodes_[6]},  // 后面
                {nodes_[0], nodes_[4], nodes_[7], nodes_[3]},  // 左面
                {nodes_[1], nodes_[2], nodes_[6], nodes_[5]}}; // 右面
    }
    
private:
    std::vector<Index> nodes_;
};

// 20节点六面体 (二次)
class Brick20 : public Volume {
public:
    Brick20(const std::array<Index,20>& nodes);
    
    ElementType type() const override { return ElementType::Brick20; }
    int num_nodes() const override { return 20; }
    int num_faces() const override { return 6; }
    int order() const override { return 2; }
    
private:
    std::vector<Index> nodes_;  // 8 corner + 12 edge
};
```

---

## 2. Mesh 类重构

```cpp
class Mesh {
public:
    Mesh() = default;
    
    // ═══ 节点管理 ═══
    Index add_node(const Vec3& coords);
    const Node& node(Index id) const;
    Node& node(Index id);
    std::size_t num_nodes() const { return nodes_.size(); }
    
    // ═══ 单元管理 ═══
    template<typename ElemT, typename... Args>
    Index add_element(Args&&... args) {
        auto elem = std::make_unique<ElemT>(std::forward<Args>(args)...);
        Index id = elements_.size();
        elements_.push_back(std::move(elem));
        
        // 按维度分类
        int dim = elements_.back()->dimension();
        element_by_dim_[dim].push_back(id);
        
        return id;
    }
    
    const Element& element(Index id) const { return *elements_[id]; }
    Element& element(Index id) { return *elements_[id]; }
    
    std::size_t num_elements() const { return elements_.size(); }
    std::size_t num_volumes() const { return element_by_dim_[3].size(); }
    std::size_t num_faces() const { return element_by_dim_[2].size(); }
    std::size_t num_edges() const { return element_by_dim_[1].size(); }
    
    // 按维度访问
    const std::vector<Index>& volumes() const { return element_by_dim_[3]; }
    const std::vector<Index>& faces() const { return element_by_dim_[2]; }
    const std::vector<Index>& edges() const { return element_by_dim_[1]; }
    
    // ═══ 材料管理 ═══
    void set_material(Index vol_id, int material_id);
    int material_id(Index vol_id) const;
    
    // 按材料获取体单元
    std::vector<Index> volumes_with_material(int material_id) const;
    
    // ═══ 区域/物理域管理 ═══
    void add_region(const std::string& name, const std::vector<Index>& elem_ids);
    const std::vector<Index>& region(const std::string& name) const;
    bool has_region(const std::string& name) const;
    
    // ═══ 边界管理 ═══
    void add_boundary(const std::string& name, const std::vector<Index>& face_ids);
    const std::vector<Index>& boundary(const std::string& name) const;
    bool has_boundary(const std::string& name) const;
    
    // ═══ 拓扑查询 ═══
    // 查找节点相邻的所有单元
    std::vector<Index> elements_containing_node(Index node_id) const;
    
    // 查找面的相邻体单元
    std::vector<Index> volumes_adjacent_to_face(Index face_id) const;
    
    // 提取体单元的边界面
    std::vector<Index> extract_boundary_faces() const;
    
    // ═══ 辅助功能 ═══
    void print_info() const;
    void validate() const;  // 检查拓扑一致性
    
private:
    // 存储
    std::vector<Node> nodes_;
    std::vector<std::unique_ptr<Element>> elements_;
    
    // 索引
    std::array<std::vector<Index>, 4> element_by_dim_;  // [0]=nodes, [1]=edges, [2]=faces, [3]=volumes
    
    // 区域管理
    std::unordered_map<std::string, std::vector<Index>> regions_;
    std::unordered_map<std::string, std::vector<Index>> boundaries_;
    
    // 材料映射
    std::unordered_map<Index, int> volume_to_material_;
};
```

---

## 3. 材料系统

```cpp
struct Material {
    int id;
    std::string name;
    
    // 通用属性
    std::unordered_map<std::string, Real> properties;
    
    // 辅助访问
    Real get(const std::string& key, Real default_val = 0.0) const {
        auto it = properties.find(key);
        return (it != properties.end()) ? it->second : default_val;
    }
    
    void set(const std::string& key, Real value) {
        properties[key] = value;
    }
};

class MaterialLibrary {
public:
    int add_material(const std::string& name) {
        int id = materials_.size();
        materials_.push_back({id, name, {}});
        return id;
    }
    
    Material& material(int id) { return materials_[id]; }
    const Material& material(int id) const { return materials_[id]; }
    
    int find_by_name(const std::string& name) const {
        for (const auto& mat : materials_) {
            if (mat.name == name) return mat.id;
        }
        return -1;
    }
    
private:
    std::vector<Material> materials_;
};
```

---

## 4. 使用示例

### 4.1 创建 2D 三角形网格

```cpp
Mesh mesh;

// 添加节点
Index n0 = mesh.add_node({0, 0, 0});
Index n1 = mesh.add_node({1, 0, 0});
Index n2 = mesh.add_node({1, 1, 0});
Index n3 = mesh.add_node({0, 1, 0});

// 添加面单元 (视为2D问题)
Index f0 = mesh.add_element<Tri3>(n0, n1, n2);
Index f1 = mesh.add_element<Tri3>(n0, n2, n3);

// 设置边界
mesh.add_boundary("left", {/* 左边界面ID */});
mesh.add_boundary("bottom", {/* 底边界面ID */});

// 设置材料
mesh.set_material(f0, 0);  // 材料0
mesh.set_material(f1, 0);
```

### 4.2 创建 3D 四面体网格

```cpp
Mesh mesh;

// 节点
Index n0 = mesh.add_node({0, 0, 0});
Index n1 = mesh.add_node({1, 0, 0});
Index n2 = mesh.add_node({0, 1, 0});
Index n3 = mesh.add_node({0, 0, 1});

// 体单元
Index v0 = mesh.add_element<Tet4>(n0, n1, n2, n3);

// 材料
MaterialLibrary matlib;
int steel_id = matlib.add_material("Steel");
matlib.material(steel_id).set("E", 2.1e11);
matlib.material(steel_id).set("nu", 0.3);

mesh.set_material(v0, steel_id);

// 提取边界面
auto boundary_faces = mesh.extract_boundary_faces();
mesh.add_boundary("external", boundary_faces);
```

### 4.3 多材料区域

```cpp
// 钢材区域
std::vector<Index> steel_volumes = {0, 1, 2};
for (Index vid : steel_volumes) {
    mesh.set_material(vid, steel_id);
}
mesh.add_region("steel_part", steel_volumes);

// 铝材区域
std::vector<Index> aluminum_volumes = {3, 4, 5};
for (Index vid : aluminum_volumes) {
    mesh.set_material(vid, aluminum_id);
}
mesh.add_region("aluminum_part", aluminum_volumes);
```

---

## 5. 迁移路径

### Phase 1: 核心类实现
1. 实现 Element 基类和子类 (Node, Edge2, Tri3, Quad4, Tet4, Brick8)
2. 重构 Mesh 类
3. 实现 MaterialLibrary

### Phase 2: 兼容层
1. 保留旧的 mesh_generator (生成新结构)
2. 适配 Assembler 使用新 Mesh

### Phase 3: 高阶单元
1. 添加 Tri6, Quad8, Tet10, Brick20
2. 实现高阶形函数

### Phase 4: 清理
1. 移除旧代码
2. 更新测试和示例

---

## 6. 优势

✅ **通用性**: 支持 1D/2D/3D 混合网格  
✅ **可扩展**: 轻松添加新单元类型  
✅ **拓扑关系**: Volume-Face-Edge-Node 层次清晰  
✅ **材料管理**: 多材料/多区域支持  
✅ **高阶单元**: 架构支持任意阶次  
✅ **工业标准**: 符合现代 FEM 软件设计模式

---

## 7. 与当前架构对比

| 特性 | 当前 Phase 2 | 新架构 |
|------|-------------|--------|
| 单元类型 | Triangle2D, Quad2D | Tri3/6, Quad4/8, Tet4/10, Brick8/20 |
| 维度 | 仅 2D | 1D/2D/3D 统一 |
| 拓扑关系 | 无 | Volume→Face→Edge→Node |
| 材料系统 | 无 | MaterialLibrary + 单元关联 |
| 区域管理 | 简单标记 | Region + Boundary 系统 |
| 高阶单元 | 不支持 | 完整支持 |

---

**建议**: 先实现基础的 Tri3, Tet4, Brick8，验证架构可行性后再扩展高阶单元。
