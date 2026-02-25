# Material 本构模型模块

## 📋 概述

材料本构模型模块提供统一的接口和框架，用于实现各种材料行为：
- **弹性**：各向同性、正交各向异性、完全各向异性
- **塑性**：J2塑性、Drucker-Prager、Mohr-Coulomb
- **粘塑性**：率相关材料
- **超弹性**：橡胶、生物组织
- **损伤**：脆性/延性损伤模型

---

## 🏗️ 架构设计

### 核心类层次

```
Material (抽象基类)
├── IsotropicElastic (各向同性弹性) ✅
├── J2Plasticity (von Mises塑性) ✅
├── OrthotropicElastic (正交各向异性) 🚧
├── DruckerPrager (岩土材料) 🚧
├── ViscoplasticMaterial (粘塑性) 🚧
└── HyperelasticMaterial (超弹性) 🚧

✅ 已实现  🚧 计划中
```

### 状态变量管理

`StateVariables` 类管理材料的历史相关量：
- 塑性应变 `plastic_strain`
- 等效塑性应变 `equiv_plastic_strain`
- 损伤变量 `damage`
- 背应力 `back_stress`（运动硬化）
- 自定义标量/张量变量（扩展字段）

---

## 📖 使用指南

### 1. 实现新材料模型

继承 `Material` 类并实现必需接口：

```cpp
#include "material/material.h"

class MyMaterial : public Material {
public:
    MyMaterial(Real E, Real nu) : Material(6) {  // 6分量 = 3D
        setParameter("E", E);
        setParameter("nu", nu);
    }
    
    // 必须实现的接口
    void computeStress(const Vector& strain_inc, 
                       Vector& stress, 
                       StateVariables& state) override {
        // TODO: 应力更新算法
    }
    
    void computeTangent(const Vector& strain,
                        DenseMatrix& D_mat,
                        const StateVariables& state) override {
        // TODO: 切线刚度矩阵
    }
    
    Real strainEnergy(const Vector& strain,
                      const StateVariables& state) const override {
        // TODO: 应变能密度
        return 0.0;
    }
    
    StateVariables createState() const override {
        return StateVariables(6);  // 匹配应变尺寸
    }
    
    std::string typeName() const override {
        return "MyMaterial";
    }
};
```

### 2. 使用材料模型

```cpp
// 创建材料实例
auto material = std::make_shared<MyMaterial>(200e3, 0.3);

// 初始化状态
StateVariables state = material->createState();

// 应力更新（增量法）
Vector strain_inc(6);
strain_inc[0] = 0.001;  // ε11 增量
Vector stress(6);
material->computeStress(strain_inc, stress, state);

// 计算切线刚度
Vector total_strain(6);  // 当前总应变
DenseMatrix D_mat;
material->computeTangent(total_strain, D_mat, state);
```

### 3. 状态变量操作

```cpp
StateVariables state(6);

// 访问预定义变量
state.equiv_plastic_strain = 0.05;
state.damage = 0.1;
state.plastic_strain[0] = 0.002;

// 扩展标量变量
state.setScalar("kappa", 100.0);  // 硬化变量
Real kappa = state.getScalar("kappa");

// 扩展张量变量
Vector alpha(6, 0.0);
state.setTensor("kinematic_hardening", alpha);
Vector retrieved = state.getTensor("kinematic_hardening");

// 检查点：序列化
std::ofstream ofs("state.bin", std::ios::binary);
state.serialize(ofs);

// 恢复：反序列化
StateVariables loaded;
std::ifstream ifs("state.bin", std::ios::binary);
loaded.deserialize(ifs);
```

---

## 🔬 Voigt 记号约定

应力/应变张量用向量表示（Voigt记号）：

### 3D (6分量)
```
σ = [σ11, σ22, σ33, σ12, σ23, σ13]ᵀ
ε = [ε11, ε22, ε33, γ12, γ23, γ13]ᵀ
```
**注意**：工程剪应变 `γij = 2εij`

### 2D平面应力/应变 (3分量)
```
σ = [σ11, σ22, σ12]ᵀ
ε = [ε11, ε22, γ12]ᵀ
```

### 2D轴对称 (4分量)
```
σ = [σrr, σzz, σθθ, σrz]ᵀ
ε = [εrr, εzz, εθθ, γrz]ᵀ
```

---

## ✅ 测试

运行材料模块测试：

```bash
cd build
./bin/fem_tests --gtest_filter="StateVariables*:Material*"
```

当前测试覆盖：
- ✅ StateVariables 构造、赋值、序列化
- ✅ Material 参数管理、状态创建
- ✅ Mock 材料类接口验证

---

## 📚 API 参考

### Material 基类

#### 核心接口（纯虚函数，必须实现）

| 方法 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `computeStress` | 应力更新 | 应变增量 | 应力、状态 |
| `computeTangent` | 切线刚度 | 当前应变 | 刚度矩阵 |
| `strainEnergy` | 应变能 | 应变、状态 | 标量能量 |
| `createState` | 状态工厂 | - | StateVariables |
| `typeName` | 类型名 | - | 字符串 |

#### 可选覆盖

| 方法 | 默认行为 | 何时覆盖 |
|------|----------|----------|
| `computeGeometricStiffness` | 返回零矩阵 | 几何非线性 |
| `initializeState` | 零初始化 | 特殊初值 |
| `validateParameters` | 无验证 | 参数约束 |

### StateVariables

| 成员 | 类型 | 用途 |
|------|------|------|
| `plastic_strain` | Vector | 塑性应变张量 |
| `equiv_plastic_strain` | Real | 等效塑性应变 |
| `damage` | Real | 损伤变量 [0,1] |
| `back_stress` | Vector | 背应力（运动硬化）|
| `scalar_vars` | map | 自定义标量 |
| `tensor_vars` | map | 自定义张量 |

---

## 🛠️ 开发路线图

### ✅ Phase 1: 基类框架（已完成）
- [x] Material 抽象基类
- [x] StateVariables 容器
- [x] 参数管理系统
- [x] 单元测试

### ✅ Phase 2: 弹性材料（已完成）
- [x] IsotropicElastic（各向同性）
- [x] 2D/3D支持，平面应力/应变
- [x] 单元测试 + 已知解验证
- [ ] OrthotropicElastic（正交各向异性，待开发）

### ✅ Phase 3: 塑性材料（已完成）
- [x] J2Plasticity（返回映射算法）
- [x] von Mises屈服准则
- [x] 等向硬化模型
- [x] 循环加载测试
- [ ] 运动硬化（待开发）
- [ ] 一致性切线刚度（简化版已实现）

### 📅 Phase 4: 高级材料
- [ ] Drucker-Prager（岩土）
- [ ] 粘塑性（Perzyna）
- [ ] 超弹性（Neo-Hookean）

### 📅 Phase 5: 集成与优化
- [ ] 与求解器集成
- [ ] 批量材料点计算（向量化）
- [ ] 性能基准测试

---

## 📝 参考文献

1. **Simo & Hughes** (1998). *Computational Inelasticity*. Springer.
2. **de Souza Neto et al.** (2008). *Computational Methods for Plasticity*. Wiley.
3. **Holzapfel** (2000). *Nonlinear Solid Mechanics*. Wiley.

---

## 👥 贡献

当前维护：Math Agent 🧮  
项目仓库：https://github.com/blog-guang/math.git

欢迎贡献新材料模型、优化和测试！
