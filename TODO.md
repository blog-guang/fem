# FEM 开发 TODO

## ✅ 已完成

### Phase 1: Mesh V2 架构基础
- [x] Element 类层次 (Node, Edge2, Tri3, Quad4, Tet4, Brick8)
- [x] Material 材料系统
- [x] Mesh 单材料域
- [x] Model 顶层容器
- [x] 移除旧架构 (fem::v2 → fem)
- [x] 基础测试 (27/27 通过)

## 🚧 进行中

### Phase 2: 核心功能实现

#### 2.1 网格生成器 ✅
- [x] `mesh_generator.h/cpp`
  - [x] `generate_unit_square_tri()` → 生成 fem::Mesh
  - [x] `generate_unit_square_quad()`
  - [x] `generate_unit_cube_tet()`
  - [x] `generate_unit_cube_brick()`
  - [x] `identify_boundaries_2d()` - 自动识别 2D 边界
  - [x] `identify_boundaries_3d()` - 自动识别 3D 边界
- [x] 测试: `test_mesh_generator.cpp` (11/11 通过) ✅

#### 2.2 装配系统 (Assembler) ✅
- [x] `assembler.h/cpp`
  - [x] 支持单个 Mesh 装配
  - [x] 支持 Model 的多 Mesh 装配
  - [x] dofs_per_node 支持 (标量/矢量场)
  - [x] 材料参数传递 (通过 Mesh::material())
  - [x] COO 装配 → CSR 转换
- [x] `boundary_condition.h/cpp`
  - [x] Dirichlet BC (支持多分量)
  - [x] 完全消去法 (主对角线置1)
  - [ ] Neumann BC (待实现)
- [x] 测试: `test_assembler.cpp` (6/6 通过) ✅

#### 2.3 物理模块重写 ✅
- [x] `physics/heat.h/cpp`
  - [x] `HeatConduction` 类
  - [x] `compute_element()` 单元计算接口
  - [x] 导热系数 k, 热源 Q 支持
  - [x] Tri3 形函数梯度计算
- [x] `physics/elasticity_v2.h/cpp`
  - [x] `Elasticity2D` 类
  - [x] 平面应力/应变本构
  - [x] B 矩阵 (应变-位移关系)
  - [x] D 矩阵 (本构矩阵)
  - [x] Tri3 单元刚度矩阵
- [x] 测试: `test_physics.cpp` (6/6 通过) ✅

#### 2.4 IO 系统 (部分完成)
- [x] `io/vtk_writer.h/cpp`
  - [x] 适配新 Mesh ✅
  - [x] 单 Mesh 输出 ✅
  - [x] 支持 Element 类型自动识别 (6种) ✅
  - [x] 节点标量场 (`add_point_scalar`) ✅
  - [x] 节点矢量场 (`add_point_vector`, 2D/3D) ✅
  - [x] **单元标量场 (`add_cell_scalar`)** ✅
  - [x] **单元矢量场 (`add_cell_vector`, 2D/3D)** ✅
  - [ ] 多 Mesh 输出 (Model 级别)
- [x] 测试: `test_io.cpp` **(16/16 通过)** ✅

#### 2.5 示例程序 ✅
- [x] `examples/poisson_2d_v2.cpp` (使用新 Assembler) ✅
- [x] `examples/heat_2d.cpp` (新 HeatConduction 模块) ✅
- [x] `examples/elasticity_2d.cpp` (新 Elasticity2D 模块) ✅
- [x] `examples/cantilever_beam.cpp` (Neumann BC 悬臂梁) ✅
- [x] **`examples/thermal_stress_2d.cpp` (热-结构耦合)** ✅
- [ ] `examples/multi_material_2d.cpp` (多材料示例) - 可选

## 📋 待规划

### Phase 3: 高级功能
- [ ] 高阶单元 (Tri6, Quad8, Tet10, Brick20)
- [ ] 自适应网格细化 (AMR)
- [x] **预条件器 ILU(0)** ✅ (2026-02-28)
  - 不完全 LU 分解（保持原稀疏模式）
  - 前向/后向替换求解
  - 性能：迭代次数减少 62%（208 → 79）
  - 加速比：2.63x（相对于 CG），2.37x（相对于 Jacobi）
  - 5/5 测试通过，性能对比示例完成
- [x] **预条件器 AMG (AMGCL)** ✅ (2026-02-28)
  - 集成 AMGCL 库（git submodule）
  - Smoothed aggregation 粗化策略
  - 使用 AMGCL_NO_BOOST（无 Boost 依赖）
  - 性能：100x100 网格，12 次迭代（CG: 208）
  - 加速比：17.33x（迭代），2.05x（时间）
  - 4/4 测试通过，大规模问题性能最优
- [x] **非线性求解器 (Newton-Raphson)** ✅ (2026-02-25)
- [ ] 瞬态分析 (时间积分)

### Phase 4: 扩展
- [ ] GPU 加速 (CUDA)
- [ ] MPI 并行
- [ ] Python 绑定
- [ ] 更多物理场 (流体, 电磁)

## 🎯 当前优先级

**✅ 已完成 (Phase 2.1-2.5, 2026-02-25):**
- ✓ Mesh V2 架构 (Element/Material/Mesh/Model)
- ✓ 数学库 (Vector, DenseMatrix, SparseMatrixCSR/COO, 格式转换)
- ✓ 求解器 (CG, PCG + Jacobi/SSOR/ILU/AMG 预条件器)
- ✓ `mesh_generator` (4种网格类型, 2D/3D边界识别)
- ✓ **Assembler** (多自由度场支持, Dirichlet BC, Neumann BC) ✅
- ✓ **physics/heat** (HeatConduction, Tri3单元) ✅
- ✓ **physics/elasticity_v2** (Elasticity2D, 平面应力/应变) ✅
- ✓ **io/vtk_writer** (单Mesh输出, 点数据, 单元数据) ✅
- ✓ **178/178 测试全部通过** ✅
- ✓ **9个示例程序验证通过** ✅
  - benchmark_preconditioners (AMG/ILU/SSOR/Jacobi/CG 性能对比)
  - Poisson, Heat, Elasticity, Cantilever Beam, **Thermal-Stress Coupling**
- ✓ GoogleTest submodule 集成
- ✓ 代码清理 (删除 1188+ 行旧代码)
- ✓ **Neumann 边界条件** (表面力、热流、边界积分) ✅
- ✓ **悬臂梁示例** (误差 2.85%) ✅
- ✓ **热-结构耦合示例** (顺序耦合) ✅

**Phase 2 完成！Phase 3 进行中**

**✅ Phase 3.1 完成 (2026-02-25):**
- ✓ **Newton-Raphson 非线性求解器** ✅
  - 非线性问题接口 (NonlinearProblem)
  - Newton-Raphson 迭代框架
  - 线搜索算法（Backtracking）
  - 6/6 测试通过
  - 几何非线性桁架示例

**✅ Phase 3.2 完成 (2026-02-28):**
- ✓ **ILU(0) 预条件器** ✅
  - ILUPreconditioner 类实现
  - 不完全 LU 分解（zero fill-in）
  - 前向/后向替换求解
  - 集成到 PCGSolver
  - 性能：100x100 网格，迭代次数 79（CG: 208）
  - 加速比：2.63x（相对于 CG），2.37x（相对于 Jacobi）
  - 5/5 测试通过
  - benchmark_preconditioners 性能对比示例

**✅ Phase 3.3 完成 (2026-02-28):**
- ✓ **AMG 预条件器（AMGCL）** ✅
  - 集成 AMGCL 库（git submodule）
  - AMGPreconditioner 类（PIMPL 模式）
  - Smoothed aggregation + SPAI0 relaxation
  - 无 Boost 依赖（AMGCL_NO_BOOST）
  - 性能：100x100 网格，迭代次数 12（CG: 208，ILU: 79）
  - 加速比：17.33x（迭代），2.05x（时间）
  - 4/4 测试通过
  - benchmark 更新（CG/Jacobi/SSOR/ILU/AMG 对比）
  
**✅ 代码重构 (2026-02-28):**
- ✓ **Linear Solver 移动到 math/** ✅
  - solver/ → math/ (solver, cg, pcg, bicgstab, newton_raphson)
  - 更合理的代码组织结构
  - 所有 #include 路径更新
  - 所有示例和测试验证通过
  
- ✓ **示例代码清理和重组** ✅
  - 删除 7 个调试/废弃文件
  - 重命名 11 个 test_* 文件（统一命名规范）
  - 示例数量: 27 → 18 个（减少 33%）
  - CMakeLists.txt 优化（使用宏，按功能分类）
  - 代码行数减少 ~2600 行
  - 178/178 测试全部通过
  
**当前任务 (Phase 3, 2026-02-28):**
1. **高阶单元** (Tri6, Quad8)
   - 更高精度
   - 更好的应力计算

2. **瞬态分析**
   - 时间积分（显式/隐式）
   - 动力学问题

**短期 (Phase 2.5-3.0):**
- 文档完善 (API 文档、教程)
- 性能分析与优化
- 高阶单元支持 (Tri6, Quad8)

**中长期 (Phase 3.0+):**
- 自适应网格细化 (AMR)
- 非线性求解器 (Newton-Raphson)
- 瞬态分析 (时间积分)
- GPU 加速 (CUDA)

---

**架构状态:**
- ✅ V2 架构完成，旧代码已清理
- ✅ 测试覆盖全面 (Core/Math/Solver/Mesh/Assembly/Physics)
- ✅ 示例程序验证 (Poisson/Heat/Elasticity)
- 🚧 IO 系统待扩展 (多Mesh, 单元数据)
