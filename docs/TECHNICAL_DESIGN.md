# FEM 多物理场仿真 App — 技术设计文档 v0.1

> 参考架构: deal.II, FEniCS/Firedrake, MOOSE, GetFEM  
> 技术栈: C++17, CMake, GoogleTest  
> 作者: Math Agent  
> 日期: 2026-02-05

---

## 1. 项目目标

构建一个**模块化、可扩展的多物理场有限元求解器**，支持：
- 二维 / 三维静态和瞬态问题
- 多种单元类型（三角形、四边形、四面体、六面体）
- 耦合多物理场求解（热传导、结构力学、流动等）
- 高效稀疏线性代数后端
- 标准 VTK 格式输出，便于可视化

初版 (v0.1) 聚焦：**2D 静态单物理场 → 2D 静态双场耦合**，后续逐步扩展。

---

## 2. 技术栈

| 层级 | 技术 |
|------|------|
| 语言 | C++17 |
| 构建 | CMake 3.20+ |
| 测试 | GoogleTest |
| 线性代数 | 自实现 CSR 稀疏矩阵 + 迭代求解器（CG, BiCGSTAB） |
| 预条件器 | Jacobi, ILU(0)（后续扩展） |
| I/O | VTK Legacy 格式输出 |
| 可选扩展 | OpenMP 并行, Eigen 稠密矩阵, PETSc 后端 |

---

## 3. 整体架构（分层）

```
┌─────────────────────────────────────────────┐
│  Application Layer (应用层)                   │
│    Driver / Config / Simulation Loop         │
├─────────────────────────────────────────────┤
│  Physics Layer (物理层)                       │
│    HeatConduction / Elasticity / ...         │
├─────────────────────────────────────────────┤
│  Solver Layer (求解层)                        │
│    LinearSolver / NonlinearSolver / Precond  │
├─────────────────────────────────────────────┤
│  Assembly Layer (装配层)                      │
│    BilinearForm / LinearForm / Assembler     │
│    BoundaryCondition                         │
├─────────────────────────────────────────────┤
│  Element Layer (单元层)                       │
│    ElementBase / Triangle / Quad / Tet / Hex |
│    Quadrature / BasisFunction                │
├─────────────────────────────────────────────┤
│  Mesh Layer (网格层)                          │
│    Mesh / Node / Cell / BoundaryMarker       │
│    MeshGenerator / MeshRefinement            │
├─────────────────────────────────────────────┤
│  Core Layer (核心层)                          │
│    Types / Logger / Config / SparseMatrix    │
│    Vector / Timer                            │
├─────────────────────────────────────────────┤
│  IO Layer (输出层)                            │
│    VTKWriter / MeshReader                    │
└─────────────────────────────────────────────┘
```

**设计原则:**
- 单向依赖: 上层依赖下层，下层不知上层
- 接口驱动: 核心类用虚基类定义接口，具体实现可替换
- 数据与逻辑分离: 数据结构独立，算法操作数据

---

## 4. 核心数据结构设计

> 设计原则: **扁平、无冗余、高效**。尽量用 POD struct + vector，减少间接引用和虚调用。

### 4.1 Mesh — 网格

最基本的数据容器。节点和单元的 id 就是它们在 vector 中的索引，不需要单独存储。

```
Mesh
  ├── coords: vector<Vec3>              // 节点坐标, coords[i] = 第 i 个节点
  ├── cells:  vector<Cell>              // 单元列表
  └── boundaries: map<string, vector<Index>>  // 边界名 → 节点索引集

Cell (轻量 POD)
  ├── type: ElementType                 // Tri2D, Quad2D, Tet3D, Hex3D
  └── nodes: array<Index, MAX_NODES>    // 连接表 (固定数组，无动态分配)
      num_nodes: uint8_t                // 实际使用节点数
```

> 用 `array<Index, MAX_NODES>` 而非 `vector<Index>`，避免每个单元一次堆分配。MAX_NODES=8 覆盖所有常见单元。

### 4.2 SparseMatrix — 稀疏矩阵

两种格式，职责明确，分阶段使用：

```
COOMatrix (装配阶段: 随机 += 快)
  ├── rows, cols: Index
  ├── row_idx: vector<Index>
  ├── col_idx: vector<Index>
  └── values: vector<Real>

CSRMatrix (求解阶段: 按行遍历快, matvec 高效)
  ├── rows: Index
  ├── row_ptr: vector<Index>    // 长度 rows+1, row_ptr[i]..row_ptr[i+1] 是第 i 行的范围
  ├── col_idx: vector<Index>
  └── values: vector<Real>
```

流程: **装配 → COO → 一次转换 → CSR → 求解**。转换时顺便对重复项求和（同一 (i,j) 多次 add_entry 会累加）。

### 4.3 Element — 有限元单元

唯一需要虚基类的地方，因为不同单元类型的形函数算法不同。接口尽量薄。

```
QuadPoint (POD)
  ├── xi: Vec3                  // 参考坐标
  └── weight: Real              // 积分权重

ElementBase (虚基类)
  ├── num_nodes(): uint8_t                    // 节点数 (Tri=3, Quad=4, ...)
  ├── quad_points(): span<const QuadPoint>    // 积分点 (静态数据, 无分配)
  ├── shape_functions(xi) → span<Real>        // N_i(xi)
  └── shape_gradients(xi) → span<Vec3>       // dN_i/dxi
```

> 初版 P1/Q1 单元: num_nodes = num_dofs。高阶单元后续再扩展。

### 4.4 Field — 物理场值

纯数据，不含逻辑。就是一个带名字的向量。

```
Field
  ├── name: string              // "temperature", "displacement_x", ...
  └── values: vector<Real>      // 场值, 索引对应 Mesh 节点编号
```

> 初版 P1 单元时，场值直接绑定到节点，size = num_nodes。不需要额外映射。

### 4.5 BoundaryCondition — 边界条件

独立施加，不存在 Field 里。初版用常数值，后续可升级。

```
BoundaryCondition
  ├── type: BCType              // Dirichlet / Neumann
  ├── boundary_name: string     // 对应 Mesh 中的 boundary 标记
  └── value: Real               // 常数值 (初版)
```

> 施加时: 查找 boundary_name 对应的节点集 → 修改方程组 K 和 F。

### 4.6 数据流

```
Mesh (坐标 + 连接表 + 边界)
      │
      ↓
ElementBase (形函数 + 积分点)
      │
      ↓ 遍历每个单元, 计算单元矩阵
      ↓
COOMatrix (装配)
      │
      ↓ 转换 + 重复项求和
      ↓
CSRMatrix ──→ BoundaryCondition 修改 ──→ Solver ──→ Field (结果)
```

数据结构总结：**6 个核心结构，全部 struct/class，无嵌套容器套容器，内存布局紧凑。**

---

## 5. 模块设计详细

### 5.1 Core (核心)

| 文件 | 职责 |
|------|------|
| `types.h` | 基础类型: Real, Index, Vec2/3, Vector, 枚举 |
| `logger.h/.cpp` | 分级日志 (DEBUG/INFO/WARN/ERROR) |
| `timer.h/.cpp` | 高精度计时 |
| `config.h/.cpp` | 参数配置 (从 JSON/YAML 读入) |

### 5.2 Mesh (网格)

| 文件 | 职责 |
|------|------|
| `mesh.h/.cpp` | Mesh 数据容器 |
| `mesh_generator.h/.cpp` | 结构化网格生成 (Unit Square, Unit Cube 等) |
| `mesh_refinement.h/.cpp` | 网格细化 (均匀 refinement) |

### 5.3 Element (单元)

| 文件 | 职责 |
|------|------|
| `element_base.h` | 虚基类接口 |
| `triangle2d.h/.cpp` | 2D 线性三角形 (P1) |
| `quad2d.h/.cpp` | 2D 双线性四边形 (Q1) |
| `quadrature.h/.cpp` | Gauss 积分点规则 |

### 5.4 BoundaryCondition (边界条件)

| 文件 | 职责 |
|------|------|
| `boundary_condition.h/.cpp` | BC 数据 + 施加到 K/F 的逻辑 |

> 初版 P1 单元: DOF = 节点, 无需独立 DOFHandler。高阶单元再引入。

### 5.5 Assembly (装配)

| 文件 | 职责 |
|------|------|
| `bilinear_form.h` | 虚基类: 单元刚度矩阵接口 |
| `linear_form.h` | 虚基类: 单元载荷向量接口 |
| `assembler.h/.cpp` | 遍历单元，装配全局 K + F |

### 5.6 Solver (求解)

| 文件 | 职责 |
|------|------|
| `sparse_matrix.h/.cpp` | COO/CSR 矩阵 + matvec |
| `linear_solver.h` | 虚基类接口 |
| `cg.h/.cpp` | Conjugate Gradient |
| `bicgstab.h/.cpp` | BiCGSTAB |
| `preconditioner.h` | 预条件器基类 |
| `jacobi.h/.cpp` | Jacobi 预条件器 |

### 5.7 Physics (物理场模块)

| 文件 | 职责 |
|------|------|
| `heat_conduction.h/.cpp` | 稳态热传导: -∇·(k∇T)=Q |
| `elasticity.h/.cpp` | 平面应力/应变静态弹性 (后续) |

每个物理场模块继承 `BilinearForm` + `LinearForm`，封装弱形式。

### 5.8 IO (输入输出)

| 文件 | 职责 |
|------|------|
| `vtk_writer.h/.cpp` | VTK Legacy 格式写入 |
| `mesh_io.h/.cpp` | 网格文件读写 |

---

## 6. 功能 Roadmap

### Phase 1: 基础框架 (当前)
- [x] 项目骨架 + CMake 构建
- [ ] Core 层 (types, logger, timer)
- [ ] Mesh 层 + 结构化网格生成器
- [ ] Element 层 (Triangle2D, Quad2D)
- [ ] Assembly 层 + BoundaryCondition
- [ ] Solver 层 (CG, BiCGSTAB)
- [ ] IO 层 (VTK)
- [ ] 示例: 2D Poisson / 稳态热传导

### Phase 2: 多物理场耦合
- [ ] 物理场模块化 (BilinearForm / LinearForm 抽象)
- [ ] 平面应力弹性模块
- [ ] 热-结构耦合 (单向/双向)
- [ ] 非线性求解器 (Newton-Raphson)

### Phase 3: 性能与扩展
- [ ] 预条件器 (Jacobi, ILU(0))
- [ ] OpenMP 并行装配 + 求解
- [ ] 网格细化 (AMR)
- [ ] 高阶单元 (P2, Q2) + DOFHandler (节点以外的自由度)

### Phase 4: 易用性与文档
- [ ] JSON 参数配置
- [ ] Python 绑定 (pybind11)
- [ ] API 文档 + 用户手册

---

## 7. 开发规范

- **命名**: 类 PascalCase, 函数/变量 snake_case, 私有成员尾下划线
- **头文件**: `#pragma once`, 自给自足
- **测试**: 每个模块对应 `test_<module>.cpp`, 覆盖核心逻辑和边界情况
- **Commit**: 遵循 Conventional Commits (`feat/fix/refactor/docs`)
- **不写代码直到明确需求**

---

## 8. 示例流程 (2D 稳态热传导)

```
1. 生成网格 (MeshGenerator)
      ↓
2. 构建 FunctionSpace (Mesh + Element → DOFHandler)
      ↓
3. 定义物理场 (HeatConduction → BilinearForm + LinearForm)
      ↓
4. 装配全局方程组 (Assembler: K, F)
      ↓
5. 施加边界条件 (BoundaryCondition 修改 K, F)
      ↓
6. 求解线性方程组 (CG / BiCGSTAB)
      ↓
7. 输出结果 (VTKWriter)
```

---

*本文档随开发进展持续迭代更新。*
