# FEM 库快速入门

## 概述

这是一个现代化的 C++ 有限元分析库，支持：
- ✅ 多种单元类型（2D/3D）
- ✅ 材料非线性（塑性）
- ✅ 高效求解器（PCG 预条件）
- ✅ 精确后处理（< 1% 误差）

## 快速开始

### 1. 编译

```bash
cd fem
mkdir build && cd build
cmake ..
make -j4
```

### 2. 运行测试

```bash
# 单元测试（169 个）
./bin/fem_tests

# 弹性力学
./bin/test_unified_elasticity

# 塑性材料（PCG 求解器）
./bin/test_plasticity_pcg

# 非线性求解（Newton-Raphson）
./bin/test_plasticity_nr_complete

# 纯剪切
./bin/test_pure_shear

# 双轴拉伸
./bin/test_biaxial_tension
```

## 基本用法

### 弹性分析

```cpp
#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "solver/pcg.h"

using namespace fem;

// 1. 创建模型
Model model("my_problem");
int mat_id = model.add_material("steel");
int mesh_id = model.add_mesh("structure", mat_id);
Mesh& mesh = model.mesh(mesh_id);

// 2. 生成网格
MeshGenerator::generate_unit_square_quad(10, 10, mesh);
MeshGenerator::identify_boundaries_2d(mesh);

// 3. 创建材料和物理模块
Real E = 200e3, nu = 0.3;  // MPa
IsotropicElastic material(E, nu, 2, true);  // 2D plane_stress
ElasticityUnified physics(&material, 2);

// 4. 装配系统
Assembler assembler(model, 2);
assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
    physics.compute_element(elem_id, m, Ke, Fe);
});

// 5. 边界条件
std::vector<DirichletBC> bcs = {
    {"left", 0, 0.0},    // 左端固定 x
    {"left", 1, 0.0},    // 左端固定 y
    {"right", 0, 0.01}   // 右端拉伸
};
assembler.apply_dirichlet(bcs);

// 6. 求解
PCGSolver solver("jacobi");
solver.set_tol(1e-8);

const auto& K = assembler.matrix();
const auto& F = assembler.rhs();
std::vector<Real> u(F.size(), 0.0);

auto result = solver.solve(K, F.raw(), u);

if (result.converged) {
    std::cout << "Solved in " << result.iterations << " iterations\n";
}
```

### 塑性分析

```cpp
#include "material/j2_plasticity.h"
#include "solver/newton_raphson.h"

// 1. 创建塑性材料
Real E = 200e3, nu = 0.3, sigma_y = 250.0, H = 0.0;
J2Plasticity material(E, nu, sigma_y, H, 2);
ElasticityUnified physics(&material, 2);

// 2. 非线性求解
NewtonRaphsonParams nr_params;
nr_params.max_iter = 20;
nr_params.tol = 1e-6;

NewtonRaphsonSolver nr_solver;
nr_solver.set_params(nr_params);

// 3. 增量加载
int num_steps = 10;
Real max_displacement = 0.01;

std::vector<Real> u(mesh.num_nodes() * 2, 0.0);

for (int step = 1; step <= num_steps; ++step) {
    Real u_applied = step * max_displacement / num_steps;
    
    // 创建非线性问题并求解
    PlasticityProblem problem(model, physics, 2, u_applied, length);
    auto result = nr_solver.solve(problem, u);
    
    if (!result.converged) {
        std::cerr << "NR failed at step " << step << "\n";
        break;
    }
    
    std::cout << "Step " << step << ": NR iter = " << result.iterations << "\n";
}
```

### 后处理

```cpp
#include "postprocess/post_processor.h"
#include "data/data_manager.h"

// 1. 创建后处理器
postprocess::PostProcessor post_processor(model, &material, 2);
data::DataManager data_manager;

// 2. 计算应变和应力（保存到高斯点）
post_processor.compute_strain_stress(u, data_manager, "strain_gp", "stress_gp");

// 3. 外插到节点
post_processor.extrapolate_to_nodes(data_manager, "stress_gp", "stress_node");

// 4. 计算 von Mises 应力
post_processor.compute_von_mises(data_manager, "stress_node", "von_mises");

// 5. 获取结果
auto* vm_field = data_manager.get_field<data::RealData>("von_mises");
const auto& vm_data = vm_field->data();

std::cout << "Max von Mises stress: " << *std::max_element(vm_data.begin(), vm_data.end()) << "\n";
```

## 常见问题

### Q: 如何选择求解器？

**小问题（< 1000 DOF）**：
- `CGSolver` - 无预条件 CG，简单

**中等问题（1000-10000 DOF）**：
- `PCGSolver("jacobi")` - Jacobi 预条件，推荐

**大问题（> 10000 DOF）**：
- `PCGSolver("ssor")` - SSOR 预条件，更快

### Q: 如何提取应力？

**方法 1：从位移计算（简单，适合单轴）**
```cpp
Real strain = displacement / length;
Real stress = E * strain;  // 弹性
```

**方法 2：使用 PostProcessor（通用）**
```cpp
post_processor.compute_strain_stress(u, data_manager);
auto* stress_field = data_manager.get_field<data::VectorData>("stress_gp");
```

### Q: Newton-Raphson 不收敛怎么办？

1. **减小增量步长**：更多小步代替少数大步
2. **改善初始猜测**：使用上一步结果
3. **使用线搜索**：`nr_params.line_search_alpha = 0.5`
4. **检查边界条件**：确保约束充分

### Q: 如何验证结果？

1. **与解析解对比**（推荐）
2. **网格细化研究**
3. **能量平衡检查**
4. **与商业软件对比**

## 性能提示

### 1. 使用预条件器

```cpp
// ❌ 慢
CGSolver solver;

// ✅ 快（45% 减少迭代）
PCGSolver solver("jacobi");
```

### 2. 选择合适的容差

```cpp
solver.set_tol(1e-8);  // 推荐
// 1e-6: 快但精度一般
// 1e-10: 慢但精度高
```

### 3. 增量加载

```cpp
// ❌ 一步到位（可能不收敛）
Real u_max = 0.1;

// ✅ 增量加载（稳定）
for (int i = 1; i <= 10; ++i) {
    Real u = i * u_max / 10;
    // solve...
}
```

## 测试结果

| 测试 | 网格 | 精度 | 收敛性 |
|------|------|------|--------|
| 单轴拉伸 | 20×4 | 0.105% | ✓ |
| 纯剪切 | 10×10 | 0.000% | ✓ |
| 双轴拉伸 | 8×8 | 0.000% | ✓ |
| 弹塑性 | 10×3 | 0.000% | ✓ |

## 参考

### 文档
- `README.md` - 项目概述
- `src/*/README.md` - 模块文档
- `examples/*.cpp` - 示例程序

### 论文/书籍
- Zienkiewicz & Taylor: "The Finite Element Method"
- de Souza Neto et al.: "Computational Plasticity"
- Simo & Hughes: "Computational Inelasticity"

### 支持
- GitHub Issues: https://github.com/blog-guang/fem
- Email: max.niu@my.swjtu.edu.cn

---
最后更新：2026-02-27
版本：1.0
