# Mesh 数据结构重构设计 V2

## 核心概念

**Mesh = 一个材料域**
- 每个 Mesh 代表一块单一材料的区域
- Mesh 内所有 Element 共享相同的材料属性
- 整个模型由多个 Mesh + 它们之间的 Interface 组成

---

## 1. Element 类层次

### 1.1 Element (抽象基类)

```cpp
namespace fem {

enum class ElementType {
    // 0D
    Node,
    // 1D
    Edge2, Edge3,
    // 2D
    Tri3, Tri6, Quad4, Quad8,
    // 3D
    Tet4, Tet10, Brick8, Brick20
};

class Element {
public:
    virtual ~Element() = default;
    
    // 基本属性
    virtual ElementType type() const = 0;
    virtual int dimension() const = 0;      // 0=Node, 1=Edge, 2=Face, 3=Volume
    virtual int num_nodes() const = 0;
    virtual int order() const = 0;          // 1=线性, 2=二次
    
    // 节点访问
    virtual const std::vector<Index>& nodes() const = 0;
    
    // 局部编号 (在 Mesh 内的索引)
    void set_local_id(Index id) { local_id_ = id; }
    Index local_id() const { return local_id_; }
    
    // 标签 (用于子区域、边界条件等)
    void set_tag(int tag) { tag_ = tag; }
    int tag() const { return tag_; }
    
protected:
    Index local_id_{0};  // 在所属 Mesh 内的局部编号
    int tag_{0};
};

}  // namespace fem
```

### 1.2 Node (0D)

```cpp
class Node : public Element {
public:
    Node(const Vec3& coords) : coords_(coords) {}
    
    ElementType type() const override { return ElementType::Node; }
    int dimension() const override { return 0; }
    int num_nodes() const override { return 1; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override {
        static std::vector<Index> dummy;
        return dummy;  // Node 本身不引用其他节点
    }
    
    const Vec3& coords() const { return coords_; }
    Vec3& coords() { return coords_; }
    
private:
    Vec3 coords_;
};
```

### 1.3 Edge (1D)

```cpp
class Edge : public Element {
public:
    int dimension() const override { return 1; }
};

// 2节点线性边
class Edge2 : public Edge {
public:
    Edge2(Index n0, Index n1) : nodes_{n0, n1} {}
    
    ElementType type() const override { return ElementType::Edge2; }
    int num_nodes() const override { return 2; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;
};

// 3节点二次边
class Edge3 : public Edge {
public:
    Edge3(Index n0, Index n1, Index n_mid) : nodes_{n0, n1, n_mid} {}
    
    ElementType type() const override { return ElementType::Edge3; }
    int num_nodes() const override { return 3; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;
};
```

### 1.4 Face (2D)

```cpp
class Face : public Element {
public:
    int dimension() const override { return 2; }
    
    // 边界边 (按逆时针或右手定则)
    virtual std::vector<std::array<Index,2>> boundary_edges() const = 0;
};

// 3节点三角形
class Tri3 : public Face {
public:
    Tri3(Index n0, Index n1, Index n2) : nodes_{n0, n1, n2} {}
    
    ElementType type() const override { return ElementType::Tri3; }
    int num_nodes() const override { return 3; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 6节点三角形
class Tri6 : public Face {
public:
    Tri6(const std::array<Index,6>& ns) 
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Tri6; }
    int num_nodes() const override { return 6; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        // 每条边3个节点
        return {{nodes_[0], nodes_[1]},  // 边0: n0-n3-n1
                {nodes_[1], nodes_[2]},  // 边1: n1-n4-n2
                {nodes_[2], nodes_[0]}}; // 边2: n2-n5-n0
    }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2, n3,n4,n5]
};

// 4节点四边形
class Quad4 : public Face {
public:
    Quad4(Index n0, Index n1, Index n2, Index n3) 
        : nodes_{n0, n1, n2, n3} {}
    
    ElementType type() const override { return ElementType::Quad4; }
    int num_nodes() const override { return 4; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[3]},
                {nodes_[3], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 8节点四边形
class Quad8 : public Face {
public:
    Quad8(const std::array<Index,8>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Quad8; }
    int num_nodes() const override { return 8; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2,n3, n4,n5,n6,n7]
};
```

### 1.5 Volume (3D)

```cpp
class Volume : public Element {
public:
    int dimension() const override { return 3; }
    
    // 边界面 (每个面由节点列表定义)
    virtual std::vector<std::vector<Index>> boundary_faces() const = 0;
};

// 4节点四面体
class Tet4 : public Volume {
public:
    Tet4(Index n0, Index n1, Index n2, Index n3)
        : nodes_{n0, n1, n2, n3} {}
    
    ElementType type() const override { return ElementType::Tet4; }
    int num_nodes() const override { return 4; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[2], nodes_[1]},  // 底面 (外法向向下)
                {nodes_[0], nodes_[1], nodes_[3]},  // 侧面1
                {nodes_[1], nodes_[2], nodes_[3]},  // 侧面2
                {nodes_[2], nodes_[0], nodes_[3]}}; // 侧面3
    }
    
private:
    std::vector<Index> nodes_;
};

// 10节点四面体
class Tet10 : public Volume {
public:
    Tet10(const std::array<Index,10>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Tet10; }
    int num_nodes() const override { return 10; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;  // 4 corner + 6 edge mid
};

// 8节点六面体
class Brick8 : public Volume {
public:
    Brick8(const std::array<Index,8>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Brick8; }
    int num_nodes() const override { return 8; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[3], nodes_[2], nodes_[1]},  // z=0
                {nodes_[4], nodes_[5], nodes_[6], nodes_[7]},  // z=1
                {nodes_[0], nodes_[1], nodes_[5], nodes_[4]},  // y=0
                {nodes_[2], nodes_[3], nodes_[7], nodes_[6]},  // y=1
                {nodes_[0], nodes_[4], nodes_[7], nodes_[3]},  // x=0
                {nodes_[1], nodes_[2], nodes_[6], nodes_[5]}}; // x=1
    }
    
private:
    std::vector<Index> nodes_;
};

// 20节点六面体
class Brick20 : public Volume {
public:
    Brick20(const std::array<Index,20>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Brick20; }
    int num_nodes() const override { return 20; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;  // 8 corner + 12 edge mid
};
```

---

## 2. Material 材料类

```cpp
class Material {
public:
    Material(int id, const std::string& name)
        : id_(id), name_(name) {}
    
    int id() const { return id_; }
    const std::string& name() const { return name_; }
    
    // 属性访问
    void set_property(const std::string& key, Real value) {
        properties_[key] = value;
    }
    
    Real property(const std::string& key, Real default_val = 0.0) const {
        auto it = properties_.find(key);
        return (it != properties_.end()) ? it->second : default_val;
    }
    
    bool has_property(const std::string& key) const {
        return properties_.find(key) != properties_.end();
    }
    
    // 常用材料属性快捷访问
    Real E() const { return property("E"); }              // 杨氏模量
    Real nu() const { return property("nu"); }            // 泊松比
    Real rho() const { return property("rho"); }          // 密度
    Real k() const { return property("k"); }              // 热传导系数
    Real cp() const { return property("cp"); }            // 比热容
    Real alpha() const { return property("alpha"); }      // 热膨胀系数
    
    void set_elastic(Real E, Real nu) {
        set_property("E", E);
        set_property("nu", nu);
    }
    
    void set_thermal(Real k, Real rho, Real cp) {
        set_property("k", k);
        set_property("rho", rho);
        set_property("cp", cp);
    }
    
private:
    int id_;
    std::string name_;
    std::unordered_map<std::string, Real> properties_;
};
```

---

## 3. Mesh 类 (单一材料域)

```cpp
class Mesh {
public:
    Mesh(const std::string& name, const Material* material)
        : name_(name), material_(material) {}
    
    // ═══ 基本信息 ═══
    const std::string& name() const { return name_; }
    const Material* material() const { return material_; }
    
    int dimension() const { return dimension_; }  // 网格的主维度 (2D or 3D)
    
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
    // 边界由低维单元定义 (2D问题边界是Edge, 3D问题边界是Face)
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
    // 可以对单元打标签，用于细分区域或施加不同条件
    void add_region(const std::string& name, const std::vector<Index>& elem_ids) {
        regions_[name] = elem_ids;
        // 同时给这些单元打上标签
        for (Index id : elem_ids) {
            elements_[id]->set_tag(regions_.size());
        }
    }
    
    const std::vector<Index>& region(const std::string& name) const {
        static std::vector<Index> empty;
        auto it = regions_.find(name);
        return (it != regions_.end()) ? it->second : empty;
    }
    
    // ═══ 拓扑操作 ═══
    // 提取边界 (自动从体单元或面单元提取外边界)
    std::vector<Index> extract_external_boundary() const;
    
    // 查找节点相邻的所有单元
    std::vector<Index> elements_containing_node(Index node_id) const;
    
    // ═══ 辅助功能 ═══
    void print_info() const {
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
    
    void validate() const;  // 检查拓扑一致性、节点连接等
    
private:
    std::string name_;
    const Material* material_;
    int dimension_{0};  // 2D or 3D
    
    // 存储
    std::vector<Node> nodes_;
    std::vector<std::unique_ptr<Element>> elements_;
    
    // 索引
    std::array<std::vector<Index>, 4> elements_by_dim_;  // [0]=unused, [1]=edges, [2]=faces, [3]=volumes
    
    // 边界和区域
    std::unordered_map<std::string, std::vector<Index>> boundaries_;
    std::unordered_map<std::string, std::vector<Index>> regions_;
};
```

---

## 4. Model 类 (顶层模型容器)

```cpp
class Model {
public:
    Model(const std::string& name) : name_(name) {}
    
    // ═══ 材料库 ═══
    int add_material(const std::string& name) {
        int id = materials_.size();
        materials_.emplace_back(id, name);
        return id;
    }
    
    Material& material(int id) { return materials_[id]; }
    const Material& material(int id) const { return materials_[id]; }
    
    int find_material(const std::string& name) const {
        for (const auto& mat : materials_) {
            if (mat.name() == name) return mat.id();
        }
        return -1;
    }
    
    // ═══ Mesh 管理 ═══
    int add_mesh(const std::string& mesh_name, int material_id) {
        int id = meshes_.size();
        meshes_.emplace_back(mesh_name, &materials_[material_id]);
        return id;
    }
    
    Mesh& mesh(int id) { return meshes_[id]; }
    const Mesh& mesh(int id) const { return meshes_[id]; }
    
    std::size_t num_meshes() const { return meshes_.size(); }
    
    // 按名字查找
    int find_mesh(const std::string& name) const {
        for (std::size_t i = 0; i < meshes_.size(); ++i) {
            if (meshes_[i].name() == name) return i;
        }
        return -1;
    }
    
    // ═══ Interface 管理 ═══
    // Interface 用于定义不同 Mesh 之间的接触、连接关系
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
        for (const auto& mesh : meshes_) count += mesh.num_nodes();
        return count;
    }
    
    std::size_t total_elements() const {
        std::size_t count = 0;
        for (const auto& mesh : meshes_) count += mesh.num_elements();
        return count;
    }
    
    void print_info() const {
        FEM_INFO("Model: " + name_);
        FEM_INFO("  Materials: " + std::to_string(materials_.size()));
        FEM_INFO("  Meshes: " + std::to_string(meshes_.size()));
        FEM_INFO("  Interfaces: " + std::to_string(interfaces_.size()));
        FEM_INFO("  Total nodes: " + std::to_string(total_nodes()));
        FEM_INFO("  Total elements: " + std::to_string(total_elements()));
        
        for (std::size_t i = 0; i < meshes_.size(); ++i) {
            FEM_INFO("  --- Mesh " + std::to_string(i) + " ---");
            meshes_[i].print_info();
        }
    }
    
private:
    std::string name_;
    std::vector<Material> materials_;
    std::vector<Mesh> meshes_;
    std::vector<Interface> interfaces_;
};
```

---

## 5. 使用示例

### 5.1 单一材料 2D 问题

```cpp
Model model("Simple 2D Problem");

// 1. 定义材料
int steel_id = model.add_material("Steel");
model.material(steel_id).set_elastic(2.1e11, 0.3);

// 2. 创建 Mesh
int mesh_id = model.add_mesh("steel_plate", steel_id);
Mesh& mesh = model.mesh(mesh_id);

// 3. 添加节点
Index n0 = mesh.add_node({0, 0, 0});
Index n1 = mesh.add_node({1, 0, 0});
Index n2 = mesh.add_node({1, 1, 0});
Index n3 = mesh.add_node({0, 1, 0});

// 4. 添加面单元 (2D 问题用 Face 作为主单元)
Index f0 = mesh.add_element<Tri3>(n0, n1, n2);
Index f1 = mesh.add_element<Tri3>(n0, n2, n3);

// 5. 定义边界 (需要先添加 Edge 单元)
Index e_left = mesh.add_element<Edge2>(n0, n3);
Index e_bottom = mesh.add_element<Edge2>(n0, n1);
Index e_right = mesh.add_element<Edge2>(n1, n2);
Index e_top = mesh.add_element<Edge2>(n2, n3);

mesh.add_boundary("left", {e_left});
mesh.add_boundary("bottom", {e_bottom});

// 6. 打印信息
model.print_info();
```

### 5.2 多材料 3D 问题

```cpp
Model model("Multi-material Structure");

// 1. 材料库
int steel_id = model.add_material("Steel");
model.material(steel_id).set_elastic(2.1e11, 0.3);
model.material(steel_id).set_property("rho", 7850);

int aluminum_id = model.add_material("Aluminum");
model.material(aluminum_id).set_elastic(7e10, 0.33);
model.material(aluminum_id).set_property("rho", 2700);

// 2. 钢材区域
int mesh1_id = model.add_mesh("steel_block", steel_id);
Mesh& mesh1 = model.mesh(mesh1_id);

// ... 添加钢材区域的节点和 Volume 单元 (Tet4, Brick8等) ...

// 3. 铝材区域
int mesh2_id = model.add_mesh("aluminum_block", aluminum_id);
Mesh& mesh2 = model.mesh(mesh2_id);

// ... 添加铝材区域的节点和 Volume 单元 ...

// 4. 定义接触界面 (如果两个区域接触)
model.add_interface("steel_aluminum_contact", mesh1_id, mesh2_id);
// ... 设置接触面上的节点/面对应关系 ...

model.print_info();
```

### 5.3 热-结构耦合多材料

```cpp
Model model("Thermal-Structural Coupling");

// 材料
int steel_id = model.add_material("Steel");
model.material(steel_id).set_elastic(2.1e11, 0.3);
model.material(steel_id).set_thermal(50, 7850, 460);  // k, rho, cp
model.material(steel_id).set_property("alpha", 1.2e-5);  // 热膨胀系数

int ceramic_id = model.add_material("Ceramic");
model.material(ceramic_id).set_elastic(3e11, 0.25);
model.material(ceramic_id).set_thermal(3, 3000, 800);
model.material(ceramic_id).set_property("alpha", 8e-6);

// 每个材料一个 Mesh
int steel_mesh = model.add_mesh("steel_part", steel_id);
int ceramic_mesh = model.add_mesh("ceramic_coating", ceramic_id);

// ... 构建网格 ...

// 热分析时遍历所有 Mesh
for (int i = 0; i < model.num_meshes(); ++i) {
    const Mesh& mesh = model.mesh(i);
    Real k = mesh.material()->k();
    // 用该材料的 k 进行热分析装配
}

// 结构分析时同样遍历
for (int i = 0; i < model.num_meshes(); ++i) {
    const Mesh& mesh = model.mesh(i);
    Real E = mesh.material()->E();
    Real nu = mesh.material()->nu();
    // 用该材料的 E, nu 进行结构分析装配
}
```

---

## 6. 与旧架构对比

| 特性 | 旧架构 | 新架构 V2 |
|------|--------|----------|
| 材料管理 | Mesh内通过tag或material_id | 每个Mesh = 一个材料域 |
| 多材料 | 一个Mesh + material_id映射 | 多个Mesh + Model管理 |
| 区域管理 | 在Mesh内分region | 自然分离：每个材料=一个Mesh |
| 接触界面 | 不支持 | Interface 定义 Mesh 间连接 |
| 语义清晰度 | 需要额外逻辑区分区域 | 物理含义直观：Mesh=材料域 |

---

## 7. 优势

✅ **物理含义清晰**: Mesh = 材料域，直观易懂  
✅ **多材料天然支持**: 每个材料一个 Mesh  
✅ **独立性**: 每个 Mesh 独立管理，便于并行处理  
✅ **接触/界面**: Interface 明确定义 Mesh 间关系  
✅ **可扩展**: 适合复杂多物理场、多材料耦合  
✅ **数据局部性**: 单个 Mesh 数据紧凑，cache 友好  

---

## 8. 实现路径

### Phase 1: 基础框架
- Element 类层次 (Node, Edge2, Tri3, Quad4, Tet4, Brick8)
- Material 类
- Mesh 类 (单材料域)
- Model 类 (多 Mesh 容器)

### Phase 2: 迁移旧代码
- 旧的 mesh_generator 生成新 Mesh 对象
- Assembler 适配新架构 (遍历 Model 的所有 Mesh)

### Phase 3: 高级功能
- 高阶单元 (Tri6, Quad8, Tet10, Brick20)
- Interface 实现 (接触、绑定)
- 拓扑查询优化

### Phase 4: 清理
- 移除旧代码
- 完善文档和示例

---

**这个设计符合你的想法吗？** 每个 Mesh 就是一块材料区域，Model 管理所有 Mesh。
