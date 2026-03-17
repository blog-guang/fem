# 自由度管理（DOF Management）与方程组装配深度调研

## 目录

1. [自由度管理概述](#1-自由度管理概述)
2. [DOF 编号策略](#2-dof-编号策略)
3. [稀疏矩阵模式（Sparsity Pattern）](#3-稀疏矩阵模式sparsity-pattern)
4. [并行计算中的 DOF 管理](#4-并行计算中的-dof-管理)
5. [边界条件处理](#5-边界条件处理)
6. [当前项目分析与改进建议](#6-当前项目分析与改进建议)

---

## 1. 自由度管理概述

### 1.1 什么是自由度？

**自由度（Degree of Freedom, DOF）** 是描述物理系统状态所需的独立变量数量。

| 物理场 | 自由度类型 | 每节点 DOF 数 |
|--------|------------|---------------|
| 标量场（温度、压力） | 标量 | 1 |
| 2D 向量场（位移） | u, v | 2 |
| 3D 向量场（位移） | u, v, w | 3 |
| 3D 向量场 + 标量（位移 + 压力） | u, v, w, p | 4 |

### 1.2 DOF 管理的核心任务

```
┌─────────────────────────────────────────────────────────────┐
│                     DOF 管理核心任务                         │
├─────────────────────────────────────────────────────────────┤
│  1. 编号（Numbering）     将物理DOF映射到矩阵索引           │
│  2. 映射（Mapping）       local DOF ↔ global DOF           │
│  3. 分配（Distribution） 并行计算中分配到各处理器           │
│  4. 索引（Indexing）      构建非零模式，支持快速装配        │
└─────────────────────────────────────────────────────────────┘
```

### 1.3 典型工作流程

```cpp
// 典型 FEM 求解器工作流程
void solve_fem_problem() {
    // 1. 创建网格
    Mesh mesh = create_mesh();
    
    // 2. 创建有限元空间（选择形函数）
    FiniteElement fe("Lagrange", order=2);
    
    // 3. DOF 编号 ← 核心步骤
    DofHandler dof_handler(mesh, fe);
    dof_handler.distribute_dofs();  // 分配全局编号
    
    // 4. 创建稀疏矩阵模式（基于 DOF 连接关系）
    SparsityPattern sp = dof_handler.make_sparsity_pattern();
    
    // 5. 初始化矩阵和向量
    SparseMatrix K(sp);    // 基于模式预分配
    Vector F(dof_handler.n_dofs());
    
    // 6. 装配
    assemble_system(K, F, dof_handler);
    
    // 7. 应用边界条件
    apply_boundary_conditions(K, F, dof_handler);
    
    // 8. 求解
    solve(K, F, u);
}
```

---

## 2. DOF 编号策略

### 2.1 基础编号方法

#### 2.1.1 节点-based 编号（当前项目使用）

```cpp
// 当前 FEM 项目的编号方式
Index global_dof(Index node_id, Index dof_index, Index dofs_per_node) {
    return node_id * dofs_per_node + dof_index;
}

// 示例：3D 弹性问题
// 节点 0: DOF 0,1,2 → u, v, w
// 节点 1: DOF 3,4,5 → u, v, w
// 节点 2: DOF 6,7,8 → u, v, w
```

**优点**：
- 实现简单
- 内存访问局部性好（同一节点的 DOF 相邻存储）

**缺点**：
- 无法区分不同类型的场（如位移 + 压力）
- 难以优化矩阵带宽

#### 2.1.2 场-based 编号

```cpp
// 按场分别编号
class FieldDofHandler {
    std::vector<Index> displacement_dofs_;  // 位移场 DOF
    std::vector<Index> pressure_dofs_;       // 压力场 DOF
    
    Index global_dof(FieldType field, Index local_index) {
        if (field == DISPLACEMENT) 
            return displacement_dofs_[local_index];
        else 
            return pressure_dofs_[local_index];
    }
};
```

### 2.2 带宽优化编号

#### 2.2.1 Cuthill-McKee 算法

**目的**：减少刚度矩阵带宽（bandwidth）

```
编号前：                      编号后（带宽优化）：
                            
0 ─── 1 ─── 2               0 ─── 3 ─── 6
│     │     │                │     │     │
3 ─── 4 ─── 5               1 ─── 4 ─── 7
│     │     │                │     │     │
6 ─── 7 ─── 8               2 ─── 5 ─── 8

带宽: 8                      带宽: 3
```

**实现**：

```cpp
// 简化的 Cuthill-McKee 编号
class BandwidthOptimizer {
public:
    std::vector<Index> optimize(const AdjacencyList& graph) {
        // 1. 找到度最小的节点作为起始
        Index start = find_min_degree_node();
        
        // 2. BFS 遍历，按层次分配编号
        std::queue<Index> q;
        q.push(start);
        visited[start] = true;
        
        while (!q.empty()) {
            Index node = q.front();
            q.pop();
            new_number[node] = current_index++;
            
            // 按度升序访问邻居
            auto neighbors = get_neighbors_sorted_by_degree(node);
            for (auto n : neighbors) {
                if (!visited[n]) {
                    visited[n] = true;
                    q.push(n);
                }
            }
        }
        
        return new_number;
    }
};
```

#### 2.2.2 逆向 Cuthill-McKee (RCM)

deal.II 和 libMesh 都使用 RCM 优化：

```cpp
// libMesh 中的带宽优化
void DofMap::compute_sparsity(const MeshBase& mesh) {
    // 使用 RCM 重新编号
    SparsityPattern::Builder builder;
    
    // 构建邻接图
    Graph graph(mesh.n_nodes());
    for (auto& elem : mesh.elements()) {
        for (auto n : elem->node_ids()) {
            graph.add_edge(elem->id(), n);
        }
    }
    
    // RCM 优化
    auto perm = reverse_cuthill_mckee(graph);
    
    // 应用新编号
    apply_new_numbering(perm);
}
```

### 2.3 编号策略对比

| 策略 | 优点 | 缺点 | 适用场景 |
|------|------|------|----------|
| 节点顺序编号 | 简单、快速 | 带宽大 | 简单问题 |
| RCM 优化 | 带宽小、矩阵稀疏 | 额外计算开销 | 大规模问题 |
| 场分离编号 | 便于多物理场 | 复杂 | 耦合问题 |
| 波浪编号 | 缓存友好 | 不适合所有网格 | 自适应问题 |

---

## 3. 稀疏矩阵模式（Sparsity Pattern）

### 3.1 为什么需要 Sparsity Pattern？

**核心问题**：在不知道非零元位置的情况下，如何预分配稀疏矩阵内存？

```
刚度矩阵的非零元位置取决于：
1. 单元的节点连接关系
2. 形函数的支撑集
3. 边界条件类型

直接动态添加元素会导致：
- 频繁内存分配
- 内存碎片化
- 性能下降
```

### 3.2 Pattern 构建方法

#### 3.2.1 基于单元的模式构建

```cpp
// deal.II 方式
SparsityPattern build_sparsity_pattern(const DofHandler& dof_handler) {
    SparsityPattern sp(dof_handler.n_dofs());
    
    // 遍历所有单元
    for (auto& cell : dof_handler.active_cells()) {
        // 获取单元的 DOF 索引
        auto dofs = cell.dof_indices();
        
        // 单元内所有 DOF 两两连接
        for (auto i : dofs) {
            for (auto j : dofs) {
                sp.add(i, j);  // 记录非零位置
            }
        }
    }
    
    sp.compress();  // 构建最终模式
    return sp;
}
```

#### 3.2.2 估计非零元数量

```cpp
// libMesh 方式：估计最大耦合数
size_t DofMap::max_couplings_between_dofs() {
    // 对于 Lagrange 单元：max_coupling = (p+1)^dim
    // p: 单元阶次, dim: 维度
    
    size_t max_coupling = 1;
    for (int d = 0; d < dim; d++) {
        max_coupling *= (fe_order + 1);
    }
    
    // 乘以每节点的 DOF 数
    max_coupling *= n_vars_per_node;
    
    return max_coupling * n_dofs();
}
```

### 3.3 稀疏矩阵格式

#### 3.3.1 COO (Coordinate Format)

```cpp
// 适合装配阶段
class SparseMatrixCOO {
    std::vector<size_t> row_;    // 行索引
    std::vector<size_t> col_;    // 列索引
    std::vector<Real> value_;    // 值
    
    void add(size_t i, size_t j, Real v) {
        row_.push_back(i);
        col_.push_back(j);
        value_.push_back(v);
    }
};
```

**特点**：
- 优点：添加元素 O(1)
- 缺点：求解时需要转换为 CSR/CSC

#### 3.3.2 CSR (Compressed Sparse Row)

```cpp
// 适合矩阵-向量乘法
class SparseMatrixCSR {
    std::vector<size_t> row_ptr_;  // 行指针 (n+1 个)
    std::vector<size_t> col_idx_;  // 列索引
    std::vector<Real> values_;     // 值
    
    // 第 i 行的非零元范围: [row_ptr_[i], row_ptr_[i+1])
};
```

**特点**：
- 优点：矩阵-向量乘 O(nnz)，内存紧凑
- 缺点：添加元素需要重新分配

#### 3.3.3 格式转换

```cpp
SparseMatrixCSR coo_to_csr(const SparseMatrixCOO& coo) {
    // 1. 统计每行非零元数
    std::vector<size_t> nnz_per_row(coo.rows(), 0);
    for (size_t i = 0; i < coo.nnz(); i++) {
        nnz_per_row[coo.row()[i]]++;
    }
    
    // 2. 计算行指针
    std::vector<size_t> row_ptr(coo.rows() + 1, 0);
    for (size_t i = 0; i < coo.rows(); i++) {
        row_ptr[i+1] = row_ptr[i] + nnz_per_row[i];
    }
    
    // 3. 填充列索引和值（需要排序/合并重复项）
    // ...
    
    return CSR;
}
```

### 3.4 主流求解器的 Pattern 策略

| 求解器 | Pattern 策略 | 特点 |
|--------|-------------|------|
| **deal.II** | 先 Pattern 再 Matrix | 内存安全，两步走 |
| **libMesh** | DofMap 自动构建 | 自动处理复杂耦合 |
| **FEniCS** | 自动推断 | 用户透明，JIT 编译 |
| **当前项目** | COO → CSR 转换 | 简单，但无预分配 |

---

## 4. 并行计算中的 DOF 管理

### 4.1 并行挑战

```
单核计算：                        并行计算（4 进程）：

┌─────────────────────┐           ┌───────┬───────┐
│     所有 DOFs       │           │ P0    │ P1    │
│     在本地          │           │       │       │
│                     │    →      ├───────┼───────┤
│                     │           │ P2    │ P3    │
│                     │           │       │       │
└─────────────────────┘           └───────┴───────┘

核心问题：
1. 如何划分网格/DOF 到各进程？
2. 如何处理跨进程边界的数据？
3. 如何高效通信局部矩阵？
```

### 4.2 DOF 分布策略

#### 4.2.1 主人-客人（Owner-Computed）规则

```cpp
// libMesh 并行 DOF 分配
enum class DofOwnerRole {
    OWNER,    // 负责计算和存储
    GHOST     // 需要同步（只读）
};

class ParallelDofMap {
    // 每个 DOF 有一个"主人"进程
    std::vector<int> dof_owner_;  // DOF ID → 进程ID
    
    // 每个进程知道自己的"客人"DOF
    std::vector<std::vector<Index>> ghost_dofs_;
    
    // 通信模式
    void distribute_dofs() {
        // 1. 基于网格分区确定主人
        for (auto& dof : all_dofs) {
            int owner = mesh_partition[dof.mesh_id()];
            dof_owner_[dof.id()] = owner;
        }
        
        // 2. 计算跨边界依赖
        for (auto& elem : boundary_elements) {
            for (auto dof : elem.dofs()) {
                if (dof_owner_[dof] != my_pid) {
                    ghost_dofs_[my_pid].push_back(dof);
                }
            }
        }
    }
};
```

#### 4.2.2 分布式稀疏矩阵

```
进程 0 的局部矩阵结构：

┌──────────────────┬────────────┐
│    本 地 区      │  客人区    │
│  (本地 DOF)      │ (Ghost)   │
├──────────────────┼────────────┤
│                  │           │
│  K_local         │  K_ghost  │
│                  │           │
└──────────────────┴────────────┘

求解时：K_local * x_local = b_local - K_ghost * x_ghost
```

### 4.3 并行装配

```cpp
// 并行装配示例（libMesh）
void ParallelAssembler::assemble() {
    // 1. 每个进程只处理本地单元
    for (auto& cell : my_elements) {
        auto Ke = compute_element_matrix(cell);
        auto dofs = cell.dof_indices();
        
        // 2. 添加到本地矩阵（带偏移）
        for (auto i : dofs) {
            for (auto j : dofs) {
                if (is_local_dof(i) && is_local_dof(j)) {
                    local_K[i_local(i)][j_local(j)] += Ke[local_i][local_j];
                }
            }
        }
    }
    
    // 3. 进程间通信：同步边界
    MPI_Allreduce(local_K, global_K, ...);
}
```

### 4.4 并行策略对比

| 策略 | 优点 | 缺点 | 框架 |
|------|------|------|------|
| 主人-客人 | 负载均衡 | 通信开销 | libMesh, deal.II |
| 分布式矩阵 | 本地求解 | 需要分区 | PETSc |
| 全局装配 | 简单 | 内存瓶颈 | 串行代码 |

---

## 5. 边界条件处理

### 5.1 Dirichlet 边界条件

#### 5.1.1 置1法（当前项目使用）

```cpp
// 刚度矩阵置1法
void apply_dirichlet_1(SparseMatrixCSR& K, Vector& F, 
                       Index dof, Real value) {
    // 1. 矩阵对角元置1
    K.set(dof, dof, 1.0);
    
    // 2. 矩阵行/列置0（可选）
    for (auto j : K.row_nz(dof)) {
        if (j != dof) {
            K.set(dof, j, 0.0);
            K.set(j, dof, 0.0);  // 对称处理
        }
    }
    
    // 3. 右端项设为已知值
    F[dof] = value;
}
```

**问题**：
- 破坏矩阵对称性
- 可能影响收敛性

#### 5.1.2 拉格朗日乘子法

```cpp
// 扩展系统
class LagrangeMultiplier {
    // 原系统: K * u = F
    // 扩展系统:
    // [K   G] [u  ]   [F  ]
    // [G^T 0] [λ  ] = [c  ]
    // 其中 G 是约束矩阵，λ 是拉格朗日乘子
    
    void expand_system() {
        n_constraints = count_dirichlet_bc();
        K_expanded.resize(n + n_constraints, n + n_constraints);
        F_expanded.resize(n + n_constraints);
        
        // 复制原矩阵
        K_expanded.block(0,0, n, n) = K;
        F_expanded.head(n) = F;
        
        // 添加约束
        for (int i = 0; i < n_constraints; i++) {
            K_expanded(n + i, constraint_dof[i]) = 1;
            K_expanded(constraint_dof[i], n + i) = 1;
            F_expanded[n + i] = constraint_value[i];
        }
    }
};
```

#### 5.1.3 惩罚法

```cpp
// 惩罚法
void apply_penalty(SparseMatrixCSR& K, Vector& F,
                   Index dof, Real value, Real penalty = 1e10) {
    // 修改对角元
    Real k_old = K(dof, dof);
    K.set(dof, dof, k_old + penalty);
    
    // 修改右端项
    F[dof] += penalty * value;
}
```

### 5.2 Neumann 边界条件

```cpp
// 自然边界条件装配
void apply_neumann(SparseMatrixCSR& K, Vector& F,
                   const Mesh& mesh, const BoundaryCondition& bc) {
    // 遍历边界单元面
    for (auto& face : boundary_faces[bc.boundary_id()]) {
        // 数值积分计算 ∫_Γ N_i * t dΓ
        auto Fe = integrate_boundary_flux(face, bc.value());
        
        // 装配到全局向量
        auto dofs = face.dof_indices();
        for (size_t i = 0; i < dofs.size(); i++) {
            F[dofs[i]] += Fe[i];
        }
    }
}
```

### 5.3 边界条件处理策略对比

| 方法 | 优点 | 缺点 | 适用场景 |
|------|------|------|----------|
| 置1法 | 简单、快速 | 破坏对称性 | 小变形、线性 |
| 惩罚法 | 保持对称 | 条件数差 | 近似约束 |
| 拉格朗日乘子 | 精确 | 增加 DOF | 复杂约束 |
| 投影法 | 无需修改矩阵 | 近似 | 大规模问题 |

---

## 6. 当前项目分析与改进建议

### 6.1 当前实现分析

```cpp
// 当前 Assembler 初始化 (assembler.cpp)
void Assembler::assemble(ElementMatrixFunc elem_func) {
    // 1. COO 格式初始化（无预分配）
    K_coo_ = SparseMatrixCOO(n_dofs_, n_dofs_);
    F_ = Vector(n_dofs_, 0.0);
    
    // 2. 单元循环装配
    for (elem : elements) {
        // 计算 Ke, Fe
        elem_func(elem_id, mesh, Ke, Fe);
        
        // 添加到 COO（动态增长）
        for (i, j : elem_dofs) {
            K_coo_.add(global_i, global_j, Ke[local_i][local_j]);
        }
    }
    
    // 3. 转换为 CSR（需要排序/合并）
    K_csr = coo_to_csr(K_coo_);
}
```

### 6.2 问题诊断

| 问题 | 影响 | 优先级 |
|------|------|--------|
| 无稀疏模式预分配 | 动态扩容，性能差 | **高** |
| DOF 编号无优化 | 矩阵带宽大 | **高** |
| 逐个应用 Dirichlet | 多次矩阵操作 | 中 |
| 无并行支持 | 无法扩展到大规模 | 低 |

### 6.3 改进方案

#### 阶段 1：Sparse Pattern（高优先级）

```cpp
class DofHandler {
public:
    // 1. DOF 编号
    void distribute_dofs(const Mesh& mesh, const FiniteElement& fe);
    
    // 2. 构建稀疏模式（两遍扫描）
    SparsityPattern make_sparsity_pattern() const {
        SparsityPattern sp(n_dofs_);
        
        // 第一遍：记录所有非零位置
        for (auto& elem : elements) {
            auto dofs = elem.dof_indices();
            for (auto i : dofs) {
                for (auto j : dofs) {
                    sp.add(i, j);
                }
            }
        }
        
        sp.compress();  // 构建 CSR 结构
        return sp;
    }
    
    // 3. 初始化矩阵
    SparseMatrix create_matrix() const {
        return SparseMatrix(pattern_);
    }
};
```

#### 阶段 2：带宽优化

```cpp
void DofHandler::optimize_numbering() {
    // 构建邻接图
    Graph g(n_nodes());
    for (auto& elem : elements) {
        for (auto ni : elem.node_ids()) {
            for (auto nj : elem.node_ids()) {
                if (ni != nj) g.add_edge(ni, nj);
            }
        }
    }
    
    // RCM 优化
    auto perm = reverse_cuthill_mckee(g);
    
    // 重新映射 DOF
    apply_permutation(perm);
}
```

#### 阶段 3：批量边界条件

```cpp
void Assembler::apply_dirichlet_batch(
    const std::vector<DirichletBC>& bcs,
    SparseMatrixCSR& K, 
    Vector& F) 
{
    // 收集所有需要修改的行
    std::vector<Index> modified_rows;
    std::vector<Real> values;
    
    for (auto& bc : bcs) {
        Index dof = bc.global_dof;
        modified_rows.push_back(dof);
        
        // 置零行（跳过对角元）
        for (auto nz : K.row_nz(dof)) {
            if (nz != dof) {
                K.set(dof, nz, 0.0);
                K.set(nz, dof, 0.0);
            }
        }
        
        // 置1对角元
        K.set(dof, dof, 1.0);
        
        // 修改右端项
        F[dof] = bc.value;
    }
}
```

### 6.4 推荐代码结构

```cpp
/**
 * 推荐的新架构
 */
class LinearSystem {
public:
    // ─── 初始化 ───
    void initialize(const Model& model, const FiniteElement& fe) {
        // 1. DOF 管理
        dof_handler_.distribute_dofs(model, fe);
        
        // 2. 构建稀疏模式
        auto pattern = dof_handler_.make_sparsity_pattern();
        
        // 3. 初始化矩阵和向量（预分配内存）
        matrix_ = std::make_unique<SparseMatrix>(pattern);
        rhs_ = Vector(dof_handler_.n_dofs());
        solution_ = Vector(dof_handler_.n_dofs());
    }
    
    // ─── 装配 ───
    void assemble() {
        // 使用 Pattern，直接添加到正确位置
        for (auto& elem : elements) {
            auto Ke = compute_element_matrix(elem);
            auto dofs = elem.dof_indices();
            matrix_->add_block(dofs, dofs, Ke);
        }
    }
    
    // ─── 边界条件 ───
    void apply_bc(const BoundaryConditions& bcs) {
        // 批量应用
        matrix_->apply_dirichlet(bcs);
        rhs_->apply_dirichlet(bcs);
    }
    
    // ─── 求解 ───
    void solve() {
        solver_->solve(*matrix_, *rhs_, *solution_);
    }
    
private:
    DofHandler dof_handler_;           // DOF 编号管理
    std::unique_ptr<SparseMatrix> matrix_;   // 稀疏矩阵
    std::unique_ptr<Vector> rhs_;           // 右端项
    std::unique_ptr<Vector> solution_;       // 解向量
};
```

---

## 7. 参考文献

1. **deal.II**: https://dealii.org/
2. **libMesh**: https://libmesh.github.io/
3. **FEniCS**: https://fenicsproject.org/
4. **PETSc**: https://www.mcs.anl.gov/petsc/
5. "The Finite Element Method: Theory, Implementation, and Applications" - Mats G. Larson, Fredrik Bengzon
6. "Iterative Methods for Sparse Linear Systems" - Yousef Saad

---

**作者**: 皮皮虾 🦐  
**日期**: 2026-03-17
