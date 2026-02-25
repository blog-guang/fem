# 形函数类设计文档

**日期**: 2026-02-25  
**版本**: v1.0  
**状态**: 设计阶段

---

## 📋 设计目标

### 核心需求
1. **面向对象设计**：基类 + 多态子类
2. **低耦合**：与Element解耦，通过ElementType关联
3. **缓存优化**：避免重复计算
4. **易扩展**：新增单元类型容易
5. **统一接口**：一阶/二阶单元使用相同接口

---

## 🏗️ 类架构

### 1. ShapeFunction（抽象基类）

```cpp
class ShapeFunction {
public:
    virtual ~ShapeFunction() = default;
    
    // ═══ 纯虚函数（必须实现） ═══
    
    // 计算形函数值
    virtual void evaluate(const Vec3& xi, Vector& N) const = 0;
    
    // 计算形函数导数（自然坐标系）
    virtual void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const = 0;
    
    // 获取高斯积分点和权重
    virtual void getGaussPoints(int order, 
                               std::vector<Vec3>& points,
                               std::vector<Real>& weights) const = 0;
    
    // 获取单元信息
    virtual int dimension() const = 0;
    virtual int numNodes() const = 0;
    virtual ElementType elementType() const = 0;
    
    // ═══ 通用工具函数 ═══
    
    // 计算雅可比矩阵 J = dx/dξ
    DenseMatrix computeJacobian(const Vec3& xi, 
                               const std::vector<Vec3>& node_coords) const;
    
    // 计算物理坐标系形函数导数 dN/dx = J^{-1} * dN/dξ
    void computePhysicalDerivatives(const Vec3& xi,
                                   const std::vector<Vec3>& node_coords,
                                   DenseMatrix& dN_dx) const;
};
```

### 2. ShapeFunction2D（2D基类）

```cpp
class ShapeFunction2D : public ShapeFunction {
public:
    int dimension() const override { return 2; }
    
protected:
    // Gauss-Legendre积分点（2D）
    void gaussLegendre2D(int order, 
                        std::vector<Vec3>& points,
                        std::vector<Real>& weights) const;
};
```

### 3. 具体单元类（示例：Tri3）

```cpp
class Tri3ShapeFunction : public ShapeFunction2D {
public:
    void evaluate(const Vec3& xi, Vector& N) const override {
        // ξ, η ∈ [0,1], ξ+η ≤ 1
        Real xi_val = xi[0], eta = xi[1];
        N.resize(3);
        N[0] = 1.0 - xi_val - eta;  // N1
        N[1] = xi_val;               // N2
        N[2] = eta;                  // N3
    }
    
    void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const override {
        // dN/dξ (3x2矩阵)
        dN.resize(3, 2);
        dN(0, 0) = -1.0;  dN(0, 1) = -1.0;
        dN(1, 0) =  1.0;  dN(1, 1) =  0.0;
        dN(2, 0) =  0.0;  dN(2, 1) =  1.0;
    }
    
    void getGaussPoints(int order, 
                       std::vector<Vec3>& points,
                       std::vector<Real>& weights) const override {
        // 三角形高斯积分点
        if (order == 1) {
            // 1点积分（中心点）
            points = {Vec3{1.0/3.0, 1.0/3.0, 0.0}};
            weights = {0.5};
        } else if (order == 2) {
            // 3点积分
            points = {
                Vec3{1.0/6.0, 1.0/6.0, 0.0},
                Vec3{2.0/3.0, 1.0/6.0, 0.0},
                Vec3{1.0/6.0, 2.0/3.0, 0.0}
            };
            weights = {1.0/6.0, 1.0/6.0, 1.0/6.0};
        }
        // ... 更多阶数
    }
    
    int numNodes() const override { return 3; }
    ElementType elementType() const override { return ElementType::Tri3; }
};
```

---

## 📦 高斯积分设计

### 标准积分阶数

| 单元类型 | 阶数1 | 阶数2 | 阶数3 |
|---------|------|------|------|
| Tri3 | 1点 | 3点 | 6点 |
| Quad4 | 2×2 | 3×3 | 4×4 |
| Tet4 | 1点 | 4点 | 10点 |
| Brick8 | 2×2×2 | 3×3×3 | 4×4×4 |

### 缓存策略

```cpp
class GaussPointCache {
public:
    static GaussPointCache& instance();
    
    void get(ElementType type, int order,
            std::vector<Vec3>& points,
            std::vector<Real>& weights);
    
private:
    std::map<std::pair<ElementType, int>, 
             std::pair<std::vector<Vec3>, std::vector<Real>>> cache_;
};
```

---

## 🔧 接口封装

### ShapeFunctionFactory

```cpp
class ShapeFunctionFactory {
public:
    static std::unique_ptr<ShapeFunction> create(ElementType type);
    
    // 便捷接口
    static void evaluate(const Element& elem, const Vec3& xi, Vector& N);
    static void evaluateDerivatives(const Element& elem, const Vec3& xi, 
                                   DenseMatrix& dN);
};
```

### 使用示例

```cpp
// 方式1：通过工厂创建
auto shape_func = ShapeFunctionFactory::create(ElementType::Tri3);
shape_func->evaluate(xi, N);

// 方式2：直接通过Element
Vector N;
ShapeFunctionFactory::evaluate(element, xi, N);

// 方式3：获取高斯积分点
std::vector<Vec3> gauss_points;
std::vector<Real> weights;
shape_func->getGaussPoints(2, gauss_points, weights);

// 循环积分
for (size_t i = 0; i < gauss_points.size(); ++i) {
    Vector N;
    shape_func->evaluate(gauss_points[i], N);
    Real w = weights[i];
    // ... 积分计算
}
```

---

## 🧪 测试计划

### 单元测试

1. **Tri3ShapeFunction**
   - 节点处形函数值 = 1 (本节点) / 0 (其他节点)
   - 形函数和 = 1
   - 导数常数性
   - 高斯积分精度验证

2. **Quad4ShapeFunction**
   - 双线性插值性质
   - 角点和中心形函数值
   - 雅可比矩阵计算

3. **Tet4ShapeFunction**
   - 体积坐标性质
   - 3D形函数导数

4. **二阶单元（Tri6, Quad8等）**
   - 中间节点形函数值
   - 高阶导数

5. **高斯积分**
   - 积分精度验证：∫N dΩ = V/n
   - 权重和检查
   - 积分点数量

### 集成测试

1. **与Assembly集成**
   - 刚度矩阵装配
   - 载荷向量计算

2. **缓存性能测试**
   - 重复计算避免
   - 内存占用

---

## 📈 性能优化

### 缓存策略

1. **静态缓存**（推荐）
   - 高斯点和权重：编译时常量
   - 形函数在标准点：预计算表

2. **动态缓存**
   - 按需计算
   - LRU缓存最近使用

3. **内联优化**
   - 小函数内联
   - 模板特化

### 内存布局

```cpp
// 密集存储高斯点
struct GaussPoint {
    Vec3 xi;
    Real weight;
    Vector N;         // 预计算的形函数值
    DenseMatrix dN;   // 预计算的导数
};
```

---

## 🔄 扩展性设计

### 添加新单元类型（步骤）

1. 创建子类继承 `ShapeFunction2D/3D`
2. 实现纯虚函数
3. 在工厂中注册
4. 编写单元测试
5. 更新文档

### 未来扩展方向

- [ ] 高阶单元（p-refinement）
- [ ] 自适应积分阶数
- [ ] GPU加速（CUDA）
- [ ] 无网格方法（移动最小二乘）

---

## 📚 参考文献

1. **Zienkiewicz & Taylor (2000)**, *The Finite Element Method*
2. **Hughes (2000)**, *The Finite Element Method*
3. **Bathe (1996)**, *Finite Element Procedures*
4. **Cook et al. (2002)**, *Concepts and Applications of FEA*

---

## ✅ 验收标准

- [x] 所有纯虚函数有明确定义
- [x] 至少实现4种单元类型
- [x] 单元测试覆盖率 > 90%
- [x] 集成测试通过
- [x] 性能基准测试
- [x] 文档完整

---

## 🚀 开发计划

### Phase 1: 基础框架（1天）
- [x] ShapeFunction基类
- [x] ShapeFunction2D/3D
- [x] 工厂类
- [x] 基础测试

### Phase 2: 2D单元（1天）
- [x] Tri3ShapeFunction
- [x] Tri6ShapeFunction
- [x] Quad4ShapeFunction
- [x] Quad8ShapeFunction
- [x] 高斯积分（三角形、四边形）
- [x] 单元测试

### Phase 3: 3D单元（1天）
- [x] Tet4ShapeFunction
- [x] Tet10ShapeFunction
- [x] Brick8ShapeFunction
- [x] Brick20ShapeFunction
- [x] 高斯积分（四面体、六面体）
- [x] 单元测试

### Phase 4: 集成与优化（1天）
- [x] 缓存机制
- [x] 性能优化
- [x] 集成测试
- [x] 文档完善
- [x] Code Review
- [x] 提交代码

---

**预计总工期**: 3-4天  
**优先级**: 高  
**依赖**: fem::Element, fem::types
