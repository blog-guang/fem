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

#### 2.1 网格生成器
- [ ] `mesh_generator.h/cpp`
  - [ ] `generate_unit_square_tri()` → 生成 fem::Mesh
  - [ ] `generate_unit_square_quad()`
  - [ ] `generate_unit_cube_tet()`
  - [ ] `generate_unit_cube_brick()`
- [ ] 测试: `test_mesh_generator.cpp`

#### 2.2 装配系统 (Assembler)
- [ ] `assembler.h/cpp`
  - [ ] 支持单个 Mesh 装配
  - [ ] 支持 Model 的多 Mesh 装配
  - [ ] dofs_per_node 支持
  - [ ] 材料参数传递 (通过 ctx)
- [ ] `boundary_condition.h/cpp`
  - [ ] Dirichlet BC (支持分量)
  - [ ] Neumann BC
- [ ] 测试: `test_assembler.cpp`

#### 2.3 物理模块重写
- [ ] `physics/heat_conduction.h/cpp`
  - [ ] HeatMaterial
  - [ ] heat_stiffness()
  - [ ] heat_load()
- [ ] `physics/elasticity.h/cpp`
  - [ ] ElasticMaterial
  - [ ] elasticity_stiffness()
  - [ ] elasticity_load()
- [ ] 测试: `test_physics.cpp`

#### 2.4 IO 系统
- [ ] `io/vtk_writer.h/cpp`
  - [ ] 适配新 Mesh
  - [ ] 支持多 Mesh 输出
  - [ ] 支持 Element 类型自动识别
- [ ] 测试: `test_io.cpp`

#### 2.5 示例程序
- [ ] `examples/poisson_2d.cpp` (使用新架构)
- [ ] `examples/heat_conduction_2d.cpp`
- [ ] `examples/elasticity_2d.cpp`
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

**立即执行 (Phase 2.1):**
1. 实现 `mesh_generator` (生成新 Mesh)
2. 编写测试验证
3. 创建一个最简示例 (验证端到端流程)

**短期 (Phase 2.2-2.3):**
- Assembler 重写
- Physics 模块适配

**中期 (Phase 2.4-2.5):**
- IO 系统
- 完整示例

---

**说明:**
- 当前架构已清理干净，无旧代码残留
- 测试覆盖: Core + Solver + Mesh/Element/Material/Model
- 下一步：从网格生成器开始，逐步重建功能链条
