# 静态结构分析实现计划

专注完善静态结构分析，对标 ANSYS Mechanical，分 7 个 Phase 逐步实现。

---

## 🎯 总体目标

**第一个月目标**：实现 ANSYS Static Structural 50% 核心功能

---

## 📊 当前状态

### ✅ 已实现（基础版）
- 线性静态求解
- 基本单元（Tri3, Quad4, Tet4, Brick8）
- 基本材料（各向同性/正交异性弹性，J2 塑性）
- 基本边界条件（固定支撑、位移、集中力）
- 基本后处理（von Mises 应力、位移、反力）
- Newton-Raphson 非线性求解器

### ❌ 缺失的关键功能
1. **Tet10 单元**（最重要！80% 工程网格）
2. **压力载荷**（最常用的载荷类型）
3. **对称边界条件**（减少模型规模）
4. **载荷步控制**（增量求解）
5. **主应力计算**（后处理必需）
6. **PARDISO 求解器**（ANSYS 默认，性能关键）
7. **接触力学**（多零件装配）

---

## 🗓️ 实施计划

### Phase 0: 增量求解框架（1周，3月第1周）

**目标**：支持载荷步和子步细分

#### 任务列表
- [ ] 设计 `LoadStep` 数据结构
  ```cpp
  struct LoadStep {
      int id;
      std::map<Index, Real> bc_start;
      std::map<Index, Real> bc_end;
      int num_substeps;
      Real time_start, time_end;
  };
  ```

- [ ] 实现载荷递增（Ramped）
  ```cpp
  for (int step = 1; step <= num_substeps; step++) {
      Real lambda = Real(step) / num_substeps;
      bc_current = bc_start + lambda * (bc_end - bc_start);
      solve(bc_current);
  }
  ```

- [ ] 实现载荷阶跃（Stepped）
  ```cpp
  // 阶跃：直接施加终止载荷
  bc_current = bc_end;
  solve(bc_current);
  ```

- [ ] 收敛历史记录
  ```cpp
  struct ConvergenceHistory {
      std::vector<Real> residuals;
      std::vector<int> iterations;
      std::vector<Real> timesteps;
  };
  ```

#### 验证
- 简单算例：悬臂梁，10 个载荷步
- 对比一次性加载 vs 增量加载
- 检查收敛历史

---

### Phase 1: Tet10 单元（2周，3月第2-3周）

**重要性**：⭐⭐⭐⭐⭐（最高优先级！）

**为什么 Tet10 如此重要**：
- Tetgen/Gmsh 自动网格的标准单元
- Tet4 太硬（shear locking），精度差
- ANSYS 中 80% 网格使用 Tet10
- 工程实际必备

#### 任务列表

1. **形函数实现**（2天）
   - [ ] 10-node 二次形函数
   ```cpp
   // 角点（顶点）
   N1 = λ1(2λ1 - 1)
   N2 = λ2(2λ2 - 1)
   N3 = λ3(2λ3 - 1)
   N4 = λ4(2λ4 - 1)
   
   // 边中点
   N5 = 4λ1λ2
   N6 = 4λ2λ3
   N7 = 4λ1λ3
   N8 = 4λ1λ4
   N9 = 4λ2λ4
   N10 = 4λ3λ4
   ```
   
   - [ ] 形函数导数（对 ξ,η,ζ）

2. **高斯积分点**（1天）
   - [ ] 4 点高斯积分（精确到 3 次多项式）
   - [ ] 权重和坐标
   ```cpp
   // Tet10 标准 4 点积分
   Vec4 xi[] = {
       {0.585410, 0.138197, 0.138197, 0.138197},
       {0.138197, 0.585410, 0.138197, 0.138197},
       {0.138197, 0.138197, 0.585410, 0.138197},
       {0.138197, 0.138197, 0.138197, 0.585410}
   };
   Real w[] = {0.25, 0.25, 0.25, 0.25};
   ```

3. **B-bar 方法**（3天）
   - [ ] 体积锁死问题理解
   - [ ] B-bar 刚度矩阵修正
   ```cpp
   // 修正 B 矩阵的体积部分
   B_bar = B_dev + (1/3) * B_vol_averaged
   ```
   - [ ] 近不可压材料测试（ν → 0.5）

4. **单元刚度矩阵装配**（2天）
   - [ ] 集成到 `ElasticityUnified`
   - [ ] 测试：单个 Tet10 单元拉伸

5. **验证**（2天）
   - [ ] NAFEMS LE10（Z-cantilever）
     - 网格：Tet10 自动网格
     - 对比解析解：δ_tip = ...
     - 目标误差：<2%
   
   - [ ] 与 Tet4 精度对比
     - 相同网格密度
     - 精度提升应 >10 倍
   
   - [ ] ANSYS 对标
     - 相同网格导入 ANSYS
     - 应力/位移误差 <1%

#### 交付
- `shape_tet10.h/cpp` 文件
- 10+ 单元测试
- 3 个验证算例报告

---

### Phase 2: 压力载荷（1周，3月第4周）

**重要性**：⭐⭐⭐⭐⭐

**应用**：压力容器、液压、气动

#### 任务列表

1. **表面单元识别**（2天）
   - [ ] 从体单元提取表面
   ```cpp
   std::vector<Surface> extract_surfaces(const Mesh& mesh) {
       // 找到只属于一个单元的面
       std::map<Face, int> face_count;
       for (elem : mesh.elements()) {
           for (face : elem.faces()) {
               face_count[face]++;
           }
       }
       return faces_with_count_1;
   }
   ```

2. **法向量计算**（1天）
   - [ ] 单元法向（叉积）
   - [ ] 节点法向（平均）
   ```cpp
   Vec3 compute_normal(Vec3 p1, Vec3 p2, Vec3 p3) {
       Vec3 v1 = p2 - p1;
       Vec3 v2 = p3 - p1;
       Vec3 n = cross(v1, v2);
       return normalize(n);
   }
   ```

3. **压力积分**（2天）
   - [ ] 面高斯积分
   - [ ] 等效节点力
   ```cpp
   // Tri3 表面：3 点高斯积分
   for (gp : gauss_points_2d) {
       Real area = 0.5 * ||v1 × v2||;
       F_e += N^T * (p * n) * area * weight;
   }
   ```

4. **大变形跟随（Follower）**（可选）
   - [ ] 法向随变形更新
   - [ ] 非线性迭代中更新

#### 验证
- [ ] 圆筒内压（解析解）
  - σ_θ = pr/t
  - σ_r = -p (内表面)
  - 误差 <1%

- [ ] 球壳内压
  - σ = pr/(2t)
  - 均匀应力场
  - 对称性验证

- [ ] ANSYS 对标：简单压力容器

---

### Phase 3: 对称边界条件（3天，4月第1周前半）

**重要性**：⭐⭐⭐⭐

**好处**：
- 减少模型规模（1/4, 1/8）
- 节省计算时间
- 节省内存

#### 任务列表

1. **对称面识别**（1天）
   - [ ] 用户指定平面（YZ, XZ, XY）
   - [ ] 自动识别对称面节点
   ```cpp
   // 例：YZ 平面（x=0）
   for (node : mesh.nodes()) {
       if (abs(node.x) < tolerance) {
           symmetry_nodes.push_back(node.id);
       }
   }
   ```

2. **法向约束施加**（1天）
   - [ ] 法向位移 = 0
   - [ ] 切向自由
   ```cpp
   // YZ 对称面：约束 DOF_X
   for (node_id : symmetry_nodes) {
       assembler.apply_constraint(node_id, DOF_X, 0.0);
   }
   ```

3. **验证**（1天）
   - [ ] 立方体拉伸（1/8 模型 vs 全模型）
   - [ ] 圆柱内压（1/4 模型）
   - [ ] 结果完全一致（<0.1% 差异）

#### 交付
- `Assembler::apply_symmetry(plane, nodes)` 接口
- 3 个验证算例

---

### Phase 4: 主应力后处理（1周，4月第1周后半）

**重要性**：⭐⭐⭐⭐

**应用**：
- 脆性材料（最大主应力准则）
- 应力分析可视化
- 安全系数计算

#### 任务列表

1. **特征值求解**（2天）
   - [ ] 选择方案：
     - Eigen 库（推荐）
     - 手写 Cardano 公式（3D 解析解）
     - QR 算法
   
   - [ ] 封装接口
   ```cpp
   void compute_principal_stress(
       const Vector& stress,  // [σxx, σyy, σzz, τxy, τyz, τxz]
       Real& sigma1,          // 最大主应力
       Real& sigma2,          // 中间主应力
       Real& sigma3,          // 最小主应力
       Vec3& dir1,            // 主方向 1
       Vec3& dir2,            // 主方向 2
       Vec3& dir3             // 主方向 3
   );
   ```

2. **后处理集成**（2天）
   - [ ] 单元中心主应力
   - [ ] 节点主应力（外推）
   - [ ] VTK 输出

3. **安全系数**（1天）
   - [ ] von Mises 准则：`SF = σ_y / σ_vm`
   - [ ] Tresca 准则：`SF = σ_y / (σ1 - σ3)`
   - [ ] 最大主应力准则：`SF = σ_y / σ1`

#### 验证
- [ ] 单轴拉伸：σ1 = σ_xx, σ2=σ3=0
- [ ] 纯剪切：σ1 = -σ3 = τ
- [ ] 双轴拉伸：σ1=σ2=σ, σ3=0

#### 交付
- `PostProcessor::compute_principal_stress()` 方法
- 安全系数计算功能
- VTK 可视化

---

### Phase 5: Multilinear Plasticity（1-2周，4月第2周）

**重要性**：⭐⭐⭐⭐

**应用**：
- 实际材料非线性行为
- 实验数据直接输入（拉伸曲线）

#### 任务列表

1. **数据结构**（1天）
   - [ ] 应力-应变曲线表格
   ```cpp
   class MultilinearIsotropic : public Material {
       std::map<Real, Real> stress_strain_curve;  // ε_p → σ_y
       
       Real get_yield_stress(Real equiv_plastic_strain) {
           // 线性插值
           auto it_upper = curve.upper_bound(eps_p);
           auto it_lower = std::prev(it_upper);
           // 插值公式
           Real sigma = linear_interpolate(
               it_lower->first, it_lower->second,
               it_upper->first, it_upper->second,
               eps_p
           );
           return sigma;
       }
   };
   ```

2. **返回映射算法扩展**（3天）
   - [ ] 迭代更新屈服应力
   ```cpp
   while (!converged) {
       sigma_y_current = get_yield_stress(eps_p_equiv);
       // Newton 迭代求解塑性乘子
       d_lambda = solve_plastic_multiplier(sigma_y_current);
       eps_p_equiv += d_lambda;
   }
   ```

3. **切线刚度**（2天）
   - [ ] 考虑硬化曲线斜率
   ```cpp
   Real H_tangent = d_sigma_y / d_eps_p;  // 表格差分
   ```

4. **验证**（2天）
   - [ ] ANSYS 对标：塑性拉伸
     - 输入：实验拉伸曲线
     - 对比：应力-应变响应
     - 误差：<5%
   
   - [ ] NAFEMS NLE1（塑性弯曲）

#### 交付
- `MultilinearIsotropic` 材料类
- 5+ 单元测试
- ANSYS 对标报告

---

### Phase 6: PARDISO 求解器集成（1周，4月第3周）

**重要性**：⭐⭐⭐⭐⭐

**原因**：
- ANSYS 默认求解器（中小规模）
- 性能优异（比 PCG 快 2-10 倍）
- 稳定性好
- 工程应用标配

#### 任务列表

1. **Intel MKL 安装**（1天）
   - [ ] 下载 Intel oneAPI Math Kernel Library
   - [ ] CMake 配置
   ```cmake
   find_package(MKL REQUIRED)
   target_link_libraries(fem PUBLIC MKL::MKL)
   ```

2. **PARDISO 封装**（2天）
   - [ ] 创建 `PARDISOSolver` 类
   ```cpp
   class PARDISOSolver : public LinearSolver {
       void* pt[64];     // PARDISO 内部指针
       int iparm[64];    // 参数数组
       int mtype = 11;   // 实对称正定
       
       SolveResult solve(const SparseMatrixCSR& K,
                        const Vector& F,
                        Vector& x) override;
   };
   ```

3. **SolverFactory 注册**（1天）
   ```cpp
   SolverFactory::registerCreator("PARDISO", 
       [](auto&, auto&) { return new PARDISOSolver(); }
   );
   ```

4. **性能测试**（2天）
   - [ ] 10k DOF：PARDISO vs PCG
   - [ ] 100k DOF：PARDISO vs PCG+AMG
   - [ ] 记录时间、内存占用

#### 验证
- ANSYS 相同算例性能对比
- 精度验证（结果完全一致）

#### 交付
- `pardiso_solver.h/cpp` 文件
- 性能测试报告
- 用户文档

---

### Phase 7: 接触力学基础（2-3周，4月第4周 - 5月第1周）

**重要性**：⭐⭐⭐⭐⭐

**应用**：
- 多零件装配
- 螺栓连接
- 压配
- 轴承

#### 任务列表

1. **Penalty Method 实现**（1周）
   - [ ] 接触检测（间隙计算）
   ```cpp
   Real gap = (x_slave - x_master) · n_master;
   if (gap < 0) {
       // 穿透，施加接触力
       F_contact = k_n * gap * n_master;
   }
   ```
   
   - [ ] 罚参数选择
   ```cpp
   k_n = alpha * E / h;  // alpha = 1~100
   ```

2. **Frictionless Contact**（3天）
   - [ ] 仅法向接触
   - [ ] 切向自由滑动

3. **Frictional Contact (Coulomb)**（1周）
   - [ ] Stick-slip 判断
   ```cpp
   Real tau_trial = k_t * delta_u_tangent;
   Real tau_max = mu * sigma_n;
   
   if (||tau_trial|| < tau_max) {
       // Stick
       tau = tau_trial;
   } else {
       // Slip
       tau = tau_max * normalize(tau_trial);
   }
   ```

4. **验证**（3天）
   - [ ] Hertz 接触（解析解）
     - 两球接触
     - 接触半径 `a`
     - 最大压力 `p_max`
     - 误差 <5%
   
   - [ ] ANSYS 对标：
     - 螺栓预紧
     - 压配
     - 摩擦滑动

#### 交付
- 接触单元实现
- Penalty/Augmented Lagrange 求解器
- Hertz 接触验证报告
- ANSYS 对标报告

---

## 📊 时间线总结

| Phase | 功能 | 工作量 | 时间 |
|-------|------|--------|------|
| 0 | 增量求解框架 | 1周 | 3月第1周 |
| 1 | Tet10 单元 | 2周 | 3月第2-3周 |
| 2 | 压力载荷 | 1周 | 3月第4周 |
| 3 | 对称边界 | 3天 | 4月第1周前半 |
| 4 | 主应力后处理 | 1周 | 4月第1周后半 |
| 5 | Multilinear Plasticity | 1-2周 | 4月第2周 |
| 6 | PARDISO | 1周 | 4月第3周 |
| 7 | 接触力学 | 2-3周 | 4月第4周 - 5月第1周 |
| **总计** | - | **8-10周** | **2个月内** |

---

## 🎯 里程碑

### M1: 核心单元完善（3月底）
- ✅ Tet10 单元
- ✅ 压力载荷
- ✅ 对称边界
- **目标**：80% 工程网格可用

### M2: 后处理增强（4月中）
- ✅ 主应力
- ✅ 安全系数
- ✅ Multilinear plasticity
- **目标**：完整的后处理功能

### M3: 性能优化（4月底）
- ✅ PARDISO 集成
- **目标**：100k DOF < 1分钟

### M4: 接触力学（5月初）
- ✅ Frictionless + Frictional contact
- **目标**：多零件装配分析

---

## 🧪 验证策略

每个 Phase 完成后必须通过：

1. **单元测试**（代码级）
   - 覆盖率 >80%
   - 边界条件测试

2. **算例验证**（功能级）
   - 解析解对比（误差 <1%）
   - NAFEMS Benchmarks

3. **ANSYS 对标**（工程级）
   - 相同网格、材料、载荷
   - 结果误差 <5%
   - 性能可比

---

## 📈 完成后的能力

实现全部 7 个 Phase 后，静态结构分析将达到：

### 功能覆盖
- ✅ ANSYS Static Structural 核心功能：**60%**
- ✅ 线性 + 非线性几何
- ✅ 弹性 + 塑性材料
- ✅ 基本接触力学
- ✅ 完整后处理

### 单元库
- ✅ Tri3, Quad4 (2D)
- ✅ Tet4, **Tet10** (3D)
- ✅ Brick8 (3D)

### 材料模型
- ✅ 各向同性/正交异性弹性
- ✅ Bilinear/Multilinear 塑性
- ✅ 等向/随动硬化

### 载荷类型
- ✅ 位移、力、压力
- ✅ 重力、温度
- ✅ 对称边界

### 接触
- ✅ Bonded
- ✅ Frictionless
- ✅ Frictional (Coulomb)

### 求解器
- ✅ CG, PCG (Jacobi/SSOR/ILU/AMG)
- ✅ **PARDISO**
- ✅ Newton-Raphson (Full/Modified)

### 后处理
- ✅ 应力/应变（总/弹/塑）
- ✅ von Mises 应力
- ✅ **主应力**
- ✅ **安全系数**
- ✅ 位移、反力

### 性能
- ✅ 100k DOF < 1分钟
- ✅ 1M DOF < 10分钟
- ✅ 并行装配（OpenMP）

---

## 🚀 开始实施

建议从 **Phase 0** 开始，按顺序实施。每完成一个 Phase，立即：

1. 编写单元测试
2. 运行验证算例
3. 对标 ANSYS
4. 提交 Git
5. 更新文档

**下一步**：开始 Phase 0 - 增量求解框架！
