# 几何非线性实现计划

## 目标

实现几何非线性分析能力，支持：
- 大变形分析（橡胶、软组织）
- 屈曲分析（薄壳、杆件）
- 接触非线性

## 理论基础

### Total Lagrangian (TL) 公式

**参考构型**：初始构型（未变形）  
**应变度量**：Green-Lagrange 应变 **E**  
**应力度量**：2nd Piola-Kirchhoff 应力 **S**

**虚功原理**：
```
δW_int = ∫ S : δE dV₀
δW_ext = ∫ f · δu dV + ∫ t · δu dA
```

**有限元离散**：
```
位移插值：u = N û
变形梯度：F = I + ∇u = I + ∂N/∂X û
Green-Lagrange 应变：E = 1/2(F^T F - I)
```

**切线刚度矩阵**：
```
K_t = K_mat + K_geo

K_mat = ∫ B_L^T D B_L dV  (材料刚度)
K_geo = ∫ G^T S G dV      (几何刚度)
```

其中：
- `B_L` = 线性应变-位移矩阵（小变形）
- `B_NL` = 非线性应变-位移矩阵（大变形）
- `G` = 几何刚度矩阵的梯度算子

### Updated Lagrangian (UL) 公式

**参考构型**：当前构型（上一增量步的变形构型）  
**应变度量**：Almansi 应变 **e**  
**应力度量**：Cauchy 应力 **σ**

**优点**：
- 更符合物理直觉
- 应力更新更简单

**缺点**：
- 需要每步更新参考构型
- 坐标系变换复杂

## 实现策略

### Phase 2.1: 运动学工具 ✅

已完成：
- Kinematics 工具类
- 变形梯度 F
- Green-Lagrange 应变 E
- Jacobian J
- 应力转换（S ↔ σ）

### Phase 2.2: 几何非线性单元（当前）

**方案 A: Updated Lagrangian（推荐用于首次实现）**
- 增量式更新
- 每步重新计算 B 矩阵
- 应力直接为 Cauchy 应力 σ

**方案 B: Total Lagrangian（理论严格）**
- 所有计算基于初始构型
- 需要应力转换 S ↔ σ
- 实现复杂度更高

**当前选择**：方案 A（Updated Lagrangian 简化版）

### Phase 2.3: 几何刚度矩阵

实现 K_geo 计算：
```cpp
K_geo = ∫ G^T σ G dV
```

其中 G 矩阵定义为：
```
G = [∂N_i/∂x  0         0
     0        ∂N_i/∂y   0
     0        0         ∂N_i/∂z]
```

### Phase 2.4: Newton-Raphson 集成

集成到现有的 Newton-Raphson 求解器：
```cpp
while (!converged) {
    // 1. 计算残差 R = F_int - F_ext
    // 2. 计算切线刚度 K_t = K_mat + K_geo
    // 3. 求解 K_t Δu = -R
    // 4. 更新 u += Δu
    // 5. 更新应力状态
}
```

## 验证案例

### 1. Cook's Membrane（经典几何非线性问题）

**问题描述**：
- 梯形悬臂板
- 端部剪切载荷
- 几何非线性明显

**参数**：
- E = 250 N/m²
- ν = 0.49999（近不可压）
- 厚度 = 1 m
- 网格：Quad4 或 Tri3

**验证指标**：
- 自由端位移 v_tip
- 与 ANSYS/文献对比

### 2. 橡胶 O-ring 压缩

**问题描述**：
- 圆环径向压缩
- 大变形 + Neo-Hookean 材料
- 接触非线性

**参数**：
- Neo-Hookean: C₁₀, D₁
- 压缩量：50%
- 对称边界条件

### 3. 薄壳屈曲

**问题描述**：
- 受轴压的圆柱壳
- 几何非线性屈曲
- 需要几何刚度矩阵

## 参考文献

1. Bathe, K.J. (2006). *Finite Element Procedures*. Chapter 6: Nonlinear Analysis.
2. Belytschko et al. (2000). *Nonlinear Finite Elements for Continua and Structures*.
3. Wriggers, P. (2008). *Nonlinear Finite Element Methods*.
4. ANSYS Theory Reference: Chapter 2.3 - Large Deformation.

## 开发进度

- [x] Phase 2.1: 运动学工具（Kinematics）
- [ ] Phase 2.2: 几何非线性单元（Updated Lagrangian）
- [ ] Phase 2.3: 几何刚度矩阵 K_geo
- [ ] Phase 2.4: Newton-Raphson 集成
- [ ] Phase 2.5: Cook's Membrane 验证
- [ ] Phase 2.6: 橡胶 O-ring 验证

**预计工作量**：1-2 周  
**当前状态**：Phase 2.1 完成，Phase 2.2 进行中
