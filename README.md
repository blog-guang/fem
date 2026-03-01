# FEM - 高性能有限元分析软件

**目标**：对标 ANSYS Mechanical 的开源有限元力学仿真软件

---

## 🎯 项目愿景

构建一个功能全面、高性能、易扩展的有限元分析软件，涵盖 ANSYS Mechanical 的核心功能：

- ✅ 静态结构分析
- ✅ 模态分析
- ✅ 瞬态动力学
- ✅ 热-结构耦合
- ✅ 非线性分析（几何 + 材料）
- ✅ 接触力学
- ✅ 屈曲分析
- ✅ 疲劳分析

---

## 📖 文档结构

完整的开发文档位于 `docs/` 目录：

### 📘 理论基础
- [01 - 有限元基本理论](docs/01-theory/01-fem-fundamentals.md)
- [02 - 材料本构模型](docs/01-theory/02-material-models.md)
- [03 - 单元技术](docs/01-theory/03-element-technology.md)
- [04 - 非线性求解](docs/01-theory/04-nonlinear-solver.md)
- [05 - 接触力学](docs/01-theory/05-contact-mechanics.md)

### 🔧 功能模块
- [Preprocessing - 前处理](docs/02-modules/01-preprocessing.md)
- [Solver - 求解器](docs/02-modules/02-solver.md)
- [Postprocessing - 后处理](docs/02-modules/03-postprocessing.md)
- [Materials - 材料系统](docs/02-modules/04-materials.md)
- [Elements - 单元库](docs/02-modules/05-elements.md)

### 💻 实现指南
- [架构设计](docs/03-implementation/01-architecture.md)
- [数据结构](docs/03-implementation/02-data-structures.md)
- [求解器实现](docs/03-implementation/03-solver-implementation.md)
- [并行化策略](docs/03-implementation/04-parallelization.md)
- [性能优化](docs/03-implementation/05-performance.md)

### 📐 分析类型
- [Static Structural - 静态结构](docs/04-analysis-types/01-static-structural.md)
- [Modal - 模态分析](docs/04-analysis-types/02-modal.md)
- [Transient - 瞬态动力学](docs/04-analysis-types/03-transient.md)
- [Thermal - 热分析](docs/04-analysis-types/04-thermal.md)
- [Thermal-Structural - 热-结构耦合](docs/04-analysis-types/05-thermal-structural.md)
- [Buckling - 屈曲分析](docs/04-analysis-types/06-buckling.md)
- [Fatigue - 疲劳分析](docs/04-analysis-types/07-fatigue.md)

### 🧪 验证案例
- [NAFEMS Benchmarks](docs/05-validation/01-nafems-benchmarks.md)
- [ANSYS 对标测试](docs/05-validation/02-ansys-comparison.md)

### 🚀 路线图
- [开发计划](docs/06-roadmap/development-plan.md)
- [功能对比表](docs/06-roadmap/feature-comparison.md)

---

## 🏗️ 当前架构

```
fem/
├── src/
│   ├── core/           # 核心工具（logger, timer, types）
│   ├── math/           # 数学库（矩阵, 向量, 求解器）
│   ├── material/       # 材料本构模型
│   ├── shape/          # 形函数与单元技术
│   ├── mesh/           # 网格数据结构
│   ├── physics/        # 物理模块（弹性, 热传导等）
│   ├── assembly/       # 装配器
│   ├── postprocess/    # 后处理
│   └── io/             # 输入输出（VTK等）
├── tests/              # 单元测试
├── examples/           # 示例程序
└── docs/               # 开发文档
```

---

## 🚀 快速开始

### 编译

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### 运行测试

```bash
./bin/fem_tests
```

### 运行示例

```bash
./bin/poisson_2d_v2      # 2D Poisson 方程
./bin/elasticity_2d      # 2D 弹性力学
./bin/cantilever_beam    # 悬臂梁 (Neumann BC)
./bin/heat_2d            # 2D 热传导
```

---

## 📊 当前功能

### ✅ 已实现

- **单元类型**：Tri3, Quad4, Tet4, Brick8
- **材料模型**：各向同性弹性, J2 塑性 (等向/随动硬化), 正交异性弹性
- **求解器**：CG, PCG (Jacobi/SSOR/ILU/AMG), BiCGSTAB, Newton-Raphson
- **分析类型**：静态结构, 热传导, 非线性材料
- **后处理**：应力/应变提取, 反力计算, VTK 输出

### 🚧 开发中

- 模态分析
- 瞬态动力学
- 接触力学
- 屈曲分析
- 热-结构耦合

### 📅 计划中

- 疲劳分析
- 复合材料
- 壳单元/梁单元
- 自适应网格
- 多物理场耦合

---

## 🤝 贡献

欢迎贡献代码、报告 Bug、提出功能建议！

---

## 📄 许可证

MIT License

---

**开发者**: 皮皮虾 🦐  
**基于**: ANSYS Mechanical 理论与工程实践
