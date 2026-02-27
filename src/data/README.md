# FEM 数据管理系统

## 概述

场数据管理系统用于存储和管理有限元网格中的各种数据，支持节点、单元、面、边、高斯点等不同位置的数据。

## 设计架构

```
BaseData (抽象基类)
    ├── FieldData<T> (模板类)
    │   ├── RealData (标量场)
    │   ├── IntData (整型场)
    │   ├── BoolData (布尔场，特化)
    │   ├── StringData (字符串场)
    │   ├── VectorData (矢量场)
    │   └── TensorData (张量场)
    └── DataManager (管理器)
```

## 核心组件

### 1. DataLocation 枚举

定义数据存储位置：

```cpp
enum class DataLocation {
    Node,        // 节点（连续场）
    Element,     // 单元中心
    Face,        // 面
    Edge,        // 边
    GaussPoint,  // 高斯积分点
    Unknown
};
```

### 2. BaseData 基类

提供通用接口：
- `resize(new_size)` - 调整大小
- `clear()` - 清空数据
- `reset()` - 重置为默认值
- `clone()` - 克隆数据
- `name()`, `location()`, `size()` - 基本信息

### 3. FieldData<T> 模板类

类型安全的数据存储：

```cpp
template<typename T>
class FieldData : public BaseData {
public:
    void set(Index index, const T& value);
    const T& get(Index index) const;
    T& get(Index index);  // 可修改版本
    
    // 高斯点专用接口
    void set_gauss_point(Index elem_id, Index gp_idx, Index n_gp, const T& value);
    const T& get_gauss_point(Index elem_id, Index gp_idx, Index n_gp) const;
    
    // 运算符重载
    T& operator[](Index index);
};
```

**支持的类型：**
- `Real` - 标量（温度、压力、密度）
- `int` - 整型（材料ID、边界标记）
- `bool` - 布尔（激活标志）
- `std::string` - 字符串（名称）
- `Vector` - 矢量（位移、速度、力）
- `DenseMatrix` - 张量（应力、应变）

### 4. DataManager 管理器

管理多个场数据：

```cpp
class DataManager {
public:
    // 添加场
    template<typename T>
    T* add_field(const std::string& name, DataLocation location, Index size = 0);
    
    // 获取场
    template<typename T>
    T* get_field(const std::string& name);
    
    // 批量操作
    void resize_all(DataLocation location, Index new_size);
    void clear_all();
    void reset_all();
    
    // 查询
    bool has_field(const std::string& name) const;
    std::vector<std::string> get_field_names() const;
};
```

## 使用示例

### 基本使用

```cpp
// 创建标量场
RealData temperature("temperature", DataLocation::Node, n_nodes, 20.0);

// 设置和获取
temperature.set(node_id, 100.0);
Real T = temperature.get(node_id);

// 使用运算符
temperature[node_id] = 100.0;
Real T = temperature[node_id];
```

### 矢量场

```cpp
// 创建矢量场
VectorData displacement("displacement", DataLocation::Node, n_nodes);

// 设置位移
Vector u(3);
u[0] = 0.01;  // u_x
u[1] = 0.02;  // u_y
u[2] = 0.00;  // u_z
displacement.set(node_id, u);

// 获取位移
const Vector& u_node = displacement.get(node_id);
```

### 张量场（高斯点）

```cpp
Index n_gp_per_elem = 4;
Index total_gp = n_elems * n_gp_per_elem;

// 创建应力场
TensorData stress("stress", DataLocation::GaussPoint, total_gp);

// 设置应力（Voigt 记号）
DenseMatrix sigma(6, 1);
sigma(0, 0) = 100.0;  // σ_xx
sigma(1, 0) = 50.0;   // σ_yy
// ...

Index elem_id = 5;
Index gp_idx = 2;
stress.set_gauss_point(elem_id, gp_idx, n_gp_per_elem, sigma);

// 获取应力
const DenseMatrix& sigma_gp = stress.get_gauss_point(elem_id, gp_idx, n_gp_per_elem);
```

### DataManager 使用

```cpp
DataManager manager;

// 注册场
manager.add_field<RealData>("temperature", DataLocation::Node, n_nodes, 20.0);
manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes);
manager.add_field<TensorData>("stress", DataLocation::GaussPoint, n_elems * n_gp);
manager.add_field<IntData>("material_id", DataLocation::Element, n_elems, 1);

// 访问场
auto* temp = manager.get_field<RealData>("temperature");
temp->set(node_id, 100.0);

auto* disp = manager.get_field<VectorData>("displacement");
Vector u(3, 0.0);
disp->set(node_id, u);

// 批量操作（网格细化）
Index new_n_nodes = 500;
manager.resize_all(DataLocation::Node, new_n_nodes);

// 重置所有高斯点数据
manager.reset_all(DataLocation::GaussPoint);
```

## 典型应用场景

### 1. 热传导分析

```cpp
DataManager manager;

// 节点场
manager.add_field<RealData>("temperature", DataLocation::Node, n_nodes, 293.15);
manager.add_field<VectorData>("heat_flux", DataLocation::Node, n_nodes);

// 单元场
manager.add_field<RealData>("conductivity", DataLocation::Element, n_elems);
manager.add_field<RealData>("heat_source", DataLocation::Element, n_elems);
```

### 2. 弹性力学分析

```cpp
DataManager manager;

// 节点场
manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes);
manager.add_field<VectorData>("reaction_force", DataLocation::Node, n_nodes);

// 高斯点场
Index n_gp = n_elems * 4;  // 每单元 4 个高斯点
manager.add_field<TensorData>("stress", DataLocation::GaussPoint, n_gp);
manager.add_field<TensorData>("strain", DataLocation::GaussPoint, n_gp);

// 单元场
manager.add_field<RealData>("element_volume", DataLocation::Element, n_elems);
manager.add_field<IntData>("material_id", DataLocation::Element, n_elems);
```

### 3. 塑性分析

```cpp
DataManager manager;

// 节点场
manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes);

// 高斯点场（需要历史变量）
Index n_gp = n_elems * n_gp_per_elem;
manager.add_field<TensorData>("stress", DataLocation::GaussPoint, n_gp);
manager.add_field<TensorData>("elastic_strain", DataLocation::GaussPoint, n_gp);
manager.add_field<TensorData>("plastic_strain", DataLocation::GaussPoint, n_gp);
manager.add_field<RealData>("equivalent_plastic_strain", DataLocation::GaussPoint, n_gp);
manager.add_field<RealData>("yield_stress", DataLocation::GaussPoint, n_gp);
```

## 特性

✅ **类型安全** - 模板化设计，编译时类型检查  
✅ **灵活存储** - 支持 Node, Element, Face, Edge, GaussPoint  
✅ **高效访问** - O(1) 索引访问  
✅ **批量操作** - 支持批量 resize, clear, reset  
✅ **易于管理** - DataManager 统一管理所有场  
✅ **可扩展** - 支持任意用户自定义类型  

## 注意事项

### BoolData 特殊处理

由于 `std::vector<bool>` 的特殊实现，`BoolData` 被特化处理：

```cpp
// ❌ 错误：不能返回 bool& 引用
bool& value = bool_field[index];

// ✅ 正确：返回 bool 值
bool value = bool_field.get(index);
bool value = bool_field[index];
```

### 高斯点索引

高斯点数据存储为线性数组：

```cpp
Index global_gp_id = elem_id * n_gp_per_elem + gp_idx;
```

推荐使用专用接口：

```cpp
stress.set_gauss_point(elem_id, gp_idx, n_gp_per_elem, value);
const auto& value = stress.get_gauss_point(elem_id, gp_idx, n_gp_per_elem);
```

### 线程安全

当前实现**不是线程安全**的。多线程环境需要外部同步。

## 测试

运行测试：

```bash
cd build
./bin/fem_tests --gtest_filter="*Data*"
```

运行示例：

```bash
./bin/test_data_manager
```

## 扩展

### 添加自定义类型

```cpp
// 自定义类型
struct MaterialState {
    Real temperature;
    Real damage;
    DenseMatrix internal_variables;
};

// 使用
using MaterialStateData = FieldData<MaterialState>;

DataManager manager;
manager.add_field<MaterialStateData>("material_state", 
                                     DataLocation::GaussPoint, 
                                     n_gp);
```

## 参考

- `examples/test_data_manager.cpp` - 完整示例
- `tests/test_data.cpp` - 单元测试
