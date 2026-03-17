# 开源有限元求解器方程组初始化调研

## 概述

调研知名开源有限元软件的线性系统（方程组）初始化架构，参考项目包括：
- deal.II (C++ 有限元库)
- libMesh (C++ 并行有限元库)
- FEniCS (Python/C++ 有限元框架)
- CalculiX (Abaqus 风格开源求解器)

---

## 1. 方程组初始化的核心问题

### 1.1 需要初始化的组件

1. **刚度矩阵 (K)** - 稀疏矩阵
2. **载荷向量 (F)** - 右端项
3. **解向量 (u)** - 待求解的未知数
4. **边界条件** - Dirichlet / Neumann

### 1.2 初始化面临的技术挑战

| 挑战 | 描述 |
|------|------|
| 稀疏性 | 大型问题可能有百万级 DOFs，需要稀疏存储 |
| 内存预分配 | 预先分配足够的空间避免动态扩容 |
| 并行初始化 | 并行计算中如何高效初始化分布式矩阵 |
| 格式转换 | COO → CSR, CRS → CCS 等格式转换 |
| 边界条件处理 | 置1法、拉格朗日乘子法、惩罚法等 |

---

## 2. 主流求解器的初始化策略

### 2.1 deal.II

**架构特点**：
- 使用 `Trilinos` 或 `PETSc` 作为底层线性代数库
- `SparsityPattern` 先于 `SparseMatrix` 创建
- 支持 Ghost Vectors（并行）

**初始化流程**：
```cpp
// 1. 创建 DOF 编号
DoFHandler<dim> dof_handler(mesh);
dof_handler.distribute_dofs(fe);

// 2. 创建稀疏Pattern（预分配非零元位置）
SparsityPattern sparsity;
sparsity.copy_from(dof_handler.max_couplings_between_dofs(),
                   dof_handler.n_dofs(),
                   dof_handler.get_dofs());

// 3. 创建矩阵
SparseMatrix<double> system_matrix;
system_matrix.reinit(sparsity);

// 4. 创建向量
Vector<double> solution, system_rhs;
solution.reinit(dof_handler.n_dofs());
system_rhs.reinit(dof_handler.n_dofs());
```

**关键设计**：
- Pattern 和 Matrix 分离：先确定非零结构，再分配内存
- `max_couplings_between_dofs()` 估计非零元数量
- 支持动态添加非零元 (Pattern 编辑器)

---

### 2.2 libMesh

**架构特点**：
- 支持 MPI 并行
- 使用 `DofMap` 管理 DOF 编号
- 可选 PETSc/Trilinos/SLEPc 求解器

**初始化流程**：
```cpp
// 1. 网格和有限元
Mesh mesh;
FEType fe_type = factory.get("LAGRANGE", SECOND);
FE<dim, ExplicitGeometry>* fe = FE<dim, ExplicitGeometry>::build(mesh, fe_type);

// 2. DOF 编号
DofMap& dof_map = system.get_dof_map();
dof_map.distribute_dofs(mesh, fe);

// 3. 创建矩阵/向量 (通过求解器)
SparseMatrix<Number>* matrix;
dof_map.attach_matrix(*matrix);
matrix->init();

// 4. 初始解向量
solution.zero();
system.update();
```

**关键设计**：
- `DofMap::distribute_dofs()` 自动处理并行 DOF 分配
- `attach_matrix()` 将矩阵与 DOFMap 关联
- `matrix->init()` 根据 DOF 结构初始化稀疏矩阵

---

### 2.3 FEniCS (Python)

**架构特点**：
- 高级 Python 接口
- 底层使用 C++ (UFL → FFC → DOLFIN)
- 自动微分 + JIT 编译

**初始化流程**：
```python
# 1. 定义函数空间
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# 2. 定义变分问题
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

# 3. 组装 (自动处理稀疏矩阵格式)
A = assemble(a)  # SparseMatrix
b = assemble(L)  # Vector

# 4. 应用边界条件
bc = DirichletBC(V, u0, "on_boundary")
bc.apply(A, b)
```

**关键设计**：
- 抽象层次高，用户无需关心矩阵格式
- UFL (Unified Form Language) 描述 PDE
- JIT 编译单元矩阵计算代码

---

### 2.4 CalculiX (Abaqus 风格)

**架构特点**：
- Fortran 编写的批量有限元求解器
- 类似于 Abaqus 的输入格式 (.inp)
- 稀疏直接求解器 (Pardiso, SPOOLES)

**初始化流程**：
```
*STATIC 分析步骤
*BOUNDARY
  1, 1, 3, 0.0  (固定节点1的所有位移)
*CLOAD
  2, 2, 100.0   (节点2施加Y方向载荷)
*ELPRINT
*Solve
```

**关键设计**：
- 输入文件驱动的装配过程
- 按单元类型批量装配
- 双精度稀疏矩阵存储

---

## 3. 当前 FEM 项目的实现

### 3.1 当前架构

```cpp
// Assembler 初始化 (assembler.cpp)
K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);  // COO 格式
F_ = Vector(n_dofs_, 0.0);                    // 右端项
is_dirichlet_dof_.resize(n_dofs_, false);     // 边界标记
```

### 3.2 装配流程

```cpp
// 单元循环 → 局部 → 全局
for (elem : elements) {
    compute_element_matrix(Ke);    // 单元刚度
    assemble_to_global(K_coo_, Fe); // 添加到 COO
}
K_csr = coo_to_csr(K_coo_);       // 转换为 CSR
```

### 3.3 边界条件处理

```cpp
// 主对角线置1法
void apply_dirichlet_single(Index global_dof, Real value) {
    K.set_diagonal(global_dof, 1.0);
    F[global_dof] = value;
}
```

---

## 4. 改进建议

### 4.1 短期改进

| 改进项 | 描述 | 优先级 |
|--------|------|--------|
| 预分配非零元 | 根据网格拓扑估计 nnz，避免动态扩容 | 高 |
| 批量 Dirichlet | 一次性应用所有边界条件 | 高 |
| 矩阵格式优化 | 评估 CSR vs CSC vs Block CSR | 中 |

### 4.2 长期改进

| 改进项 | 描述 | 优先级 |
|--------|------|--------|
| DOF 编号优化 | 使用 Cuthill-McKee 带宽优化 | 中 |
| 并行装配 | MPI 并行单元循环 | 中 |
| 外部求解器 | 集成 PETSc/Trilinos | 低 |

### 4.3 推荐的初始化流程

```cpp
class LinearSystem {
public:
    void initialize(const Model& model, const FiniteElement& fe) {
        // 1. DOF 编号
        dof_map_.distribute_dofs(model, fe);
        
        // 2. 估计非零元数量
        size_t nnz_estimate = dof_map_.estimate_nnz();
        
        // 3. 创建稀疏矩阵（预分配）
        matrix_ = SparseMatrixCSR(dof_map_.n_dofs());
        matrix_.reserve(nnz_estimate);
        
        // 4. 创建向量
        rhs_ = Vector(dof_map_.n_dofs());
        solution_ = Vector(dof_map_.n_dofs());
        
        // 5. 标记边界条件
        boundary_manager_.prepare(*dof_map_);
    }
    
    void assemble() {
        // 遍历单元 → 计算 → 装配
    }
    
    void apply_bc() {
        boundary_manager_.apply(*matrix_, *rhs_);
    }
    
    void solve() {
        solver_->solve(matrix_, rhs_, solution_);
    }
};
```

---

## 5. 参考资料

1. deal.II 官方文档: https://dealii.org/
2. libMesh 官方文档: https://libmesh.github.io/
3. FEniCS 官方文档: https://fenicsproject.org/
4. PETSc 文档: https://www.mcs.anl.gov/petsc/
5. "Finite Element Method: Linear Static and Dynamic" - J. N. Reddy

---

**作者**: 皮皮虾 🦐  
**日期**: 2026-03-17
