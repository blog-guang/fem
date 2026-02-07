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

#### 2.5 示例程序 (部分完成)
- [x] `examples/poisson_2d_v2.cpp` (使用新 Assembler) ✅
- [x] `examples/heat_2d.cpp` (新 HeatConduction 模块) ✅
- [x] `examples/elasticity_2d.cpp` (新 Elasticity2D 模块) ✅
- [ ] `examples/multi_material_2d.cpp` (多材料示例)
- [ ] `examples/thermal_stress_2d.cpp` (热-结构耦合)

## 📋 待规划

### Phase 3: 高级功能
- [ ] 高阶单元 (Tri6, Quad8, Tet10, Brick20)
- [ ] 自适应网格细化 (AMR)
- [ ] 预条件器 (ILU0, AMG)
- [ ] 非线性求解器 (Newton-Raphson)
- [ ] 瞬态分析 (时间积分)

### Phase 4: 扩展
- [ ] GPU 加速 (CUDA)
- [ ] MPI 并行
- [ ] Python 绑定
- [ ] 更多物理场 (流体, 电磁)

## 🎯 当前优先级

**✅ 已完成 (Phase 2.1-2.5, 2026-02-07):**
- ✓ Mesh V2 架构 (Element/Material/Mesh/Model)
- ✓ 数学库 (Vector, DenseMatrix, SparseMatrixCSR/COO, 格式转换)
- ✓ 求解器 (CG + Jacobi预条件器)
- ✓ `mesh_generator` (4种网格类型, 2D/3D边界识别)
- ✓ **Assembler** (多自由度场支持, Dirichlet BC, Neumann BC) ✅
- ✓ **physics/heat** (HeatConduction, Tri3单元) ✅
- ✓ **physics/elasticity_v2** (Elasticity2D, 平面应力/应变) ✅
- ✓ **io/vtk_writer** (单Mesh输出, 点数据, 单元数据) ✅
- ✓ **86/86 测试全部通过** ✅
- ✓ 7个示例程序验证通过 ✅
- ✓ GoogleTest submodule 集成
- ✓ 代码清理 (删除 1188+ 行旧代码)
- ✓ **Neumann 边界条件** (表面力、热流、边界积分) ✅
- ✓ **悬臂梁示例** (误差 2.85%) ✅

**当前任务 (Phase 2.6+, 2026-02-07):**
1. **可选扩展**
   - [ ] 多 Mesh 输出 (Model 级别) - 留待多材料示例时实现

2. **Neumann 边界条件** (Phase 2.5)
   - [ ] 自然边界条件接口
   - [ ] 边界积分
   - [ ] 表面力/热流支持

3. **更多示例** (Phase 2.5)
   - [ ] `multi_material_2d.cpp`
   - [ ] `thermal_stress_2d.cpp`

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
