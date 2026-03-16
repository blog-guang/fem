# FEM 重构计划 - 基于 libMesh

## 背景

当前 FEM 项目是一个从零开始实现的有限元分析软件，已实现：
- 基础单元（Tri3, Quad4, Tet4, Brick8）
- 材料模型（线弹性、J2 塑性、Neo-Hookean）
- 线性求解器（CG, PCG, BiCGSTAB）
- Newton-Raphson 非线性求解器
- VTK 后处理输出

## 重构目标

基于 libMesh 库重构，利用其成熟的基础设施：
- 网格管理 (Mesh)
- DOF 编号系统
- 单元系统 (Elem)
- 求解器接口 (PETSc/Trilinos)
- 并行计算
- 自适应网格细化 (AMR)

## 重构策略

### 阶段 1: 依赖集成
- [ ] 添加 libMesh 作为子模块/依赖
- [ ] 配置 CMake 找到 libMesh
- [ ] 验证基础编译

### 阶段 2: 核心架构迁移
- [ ] 网格层：使用 libMesh::Mesh 替代自研 Mesh
- [ ] DOF 系统：使用 libMesh::DofMap
- [ ] 单元库：复用 libMesh 内置单元

### 阶段 3: 物理模块重构
- [ ] 弹性力学：使用 libMesh::ElasticitySystem
- [ ] 热传导：使用 libMesh::HeatSystem
- [ ] 自研本构模型：通过 libMesh::Material 集成

### 阶段 4: 求解器迁移
- [ ] 线性求解：集成 PETSc/Trilinos
- [ ] 非线性求解：使用 libMesh::NonlinearImplicitSystem
- [ ] 特征值求解：使用 SLEPc 接口

### 阶段 5: 高级功能
- [ ] 并行计算
- [ ] 自适应网格细化 (AMR)
- [ ] 多物理场耦合

## 目录结构重构

```
fem/
├── src/
│   ├── core/           # 保留：日志、计时器、基础类型
│   ├── material/       # 保留：本构模型（集成到 libMesh）
│   ├── physics/       # 保留：物理模块（封装 libMesh 系统）
│   ├── assembly/      # 保留：装配逻辑（适配 libMesh）
│   ├── postprocess/   # 保留：后处理（封装 libMesh::ExodusII）
│   ├── io/            # 保留：I/O（使用 libMesh::IO）
│   ├── solver/        # 保留：高级求解算法
│   └── examples/      # 示例程序
├── third_party/
│   └── libmesh/      # libMesh 子模块
├── docs/              # 文档
└── tests/             # 测试
```

## 关键设计决策

### 1. 混合架构
不完全替换自研代码，而是封装 libMesh：
- 内部使用 libMesh 的基础设施
- 对外保持一致的 API 接口
- 保留核心算法实现（本构模型、求解策略）

### 2. 本构模型集成
libMesh 的材料系统相对简单，保留自研的高级本构模型：
- J2 塑性（等向/随动硬化）
- Neo-Hookean 超弹性
- 通过 libMesh::Material 基类集成

### 3. 求解器选择
- 默认使用 PETSc（推荐）
- 支持 Trilinos
- 保留自研求解器作为备选

## 实施步骤

### Step 1: 环境准备
```bash
# 克隆 libMesh
git submodule add https://github.com/libMesh/libmesh.git third_party/libmesh
cd third_party/libmesh
git submodule update --init --recursive

# 配置（需要 PETSc）
../configure --prefix=/usr/local --with-petsc-dir=$PETSC_DIR
make -j$(nproc)
```

### Step 2: 最小可行产品
1. 创建封装 libMesh::Mesh 的 Wrapper
2. 配置 CMake 链接 libMesh
3. 编译通过后运行简单示例

### Step 3: 渐进迁移
按模块逐步迁移，每迁移一个模块后确保测试通过。

## 风险与挑战

1. **编译复杂度**: libMesh 依赖较多，配置复杂
2. **API 学习曲线**: libMesh API 庞大，需要深入理解
3. **版本兼容性**: 确保 libMesh 版本兼容

## 参考资料

- libMesh 官方文档: https://libmesh.github.io/
- libMesh GitHub: https://github.com/libMesh/libmesh
- PETSc: https://www.mcs.anl.gov/petsc/

---

**作者**: 皮皮虾 🦐  
**日期**: 2026-03-16
