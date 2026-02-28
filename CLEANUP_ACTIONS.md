# FEM 库清理操作记录

## 待删除的文件（调试/废弃）

### Examples 中的调试文件
1. `examples/elasticity_debug.cpp` - 调试代码，功能已被 elasticity_2d.cpp 覆盖
2. `examples/elasticity_full_debug.cpp` - 调试代码，功能已被 elasticity_2d.cpp 覆盖
3. `examples/elasticity_simple.cpp` - 简化版，功能已被 elasticity_2d.cpp 覆盖
4. `examples/cantilever_material.cpp` - 已注释掉，功能由 cantilever_beam.cpp 提供
5. `examples/test_stress_extraction_debug.cpp` - 调试代码

## 需要重新分类的文件（test_* 系列）

这些文件以 `test_` 开头，但实际上是功能验证示例，应保留在 examples/ 但重命名：

### 保留并重命名的示例
1. `test_unified_heat.cpp` → 保留（演示统一物理模块）
2. `test_unified_elasticity.cpp` → 保留（演示统一物理模块）
3. `test_elasticity_with_plasticity.cpp` → `plasticity_demo.cpp`（演示塑性材料）
4. `test_data_manager.cpp` → `data_manager_demo.cpp`（演示数据管理）
5. `test_plasticity_verification.cpp` → `plasticity_verification.cpp`（验证塑性模型）
6. `test_plasticity_pcg.cpp` → `plasticity_pcg_demo.cpp`（演示 PCG 求解塑性）
7. `test_plasticity_nr_complete.cpp` → `plasticity_nr_demo.cpp`（演示 NR 求解塑性）
8. `test_pure_shear.cpp` → `pure_shear_demo.cpp`（纯剪切测试）
9. `test_biaxial_tension.cpp` → `biaxial_tension_demo.cpp`（双轴拉伸测试）
10. `test_incremental_postprocess.cpp` → `incremental_postprocess_demo.cpp`（增量后处理）
11. `test_reaction_forces.cpp` → `reaction_forces_demo.cpp`（反力提取）

### 可以删除的重复测试
1. `test_plasticity_simple.cpp` - 与 plasticity_pcg 重复
2. `test_plasticity_bending.cpp` - 不完整，可删除

## 核心示例（保留）

1. ✅ `poisson_2d_v2.cpp` - Poisson 方程
2. ✅ `heat_2d.cpp` - 热传导
3. ✅ `elasticity_2d.cpp` - 线弹性
4. ✅ `cantilever_beam.cpp` - 悬臂梁（Neumann BC）
5. ✅ `thermal_stress_2d.cpp` - 热-结构耦合
6. ✅ `nonlinear_truss.cpp` - 非线性求解
7. ✅ `benchmark_preconditioners.cpp` - 预条件器对比
8. ✅ `material_demo.cpp` - 材料本构演示
9. ✅ `material_simple_test.cpp` - 材料简单测试

## 执行计划

### Step 1: 删除调试文件
```bash
rm examples/elasticity_debug.cpp
rm examples/elasticity_full_debug.cpp
rm examples/elasticity_simple.cpp
rm examples/cantilever_material.cpp
rm examples/test_stress_extraction_debug.cpp
rm examples/test_plasticity_simple.cpp
rm examples/test_plasticity_bending.cpp
```

### Step 2: 重命名 test_* 文件
```bash
mv examples/test_elasticity_with_plasticity.cpp examples/plasticity_demo.cpp
mv examples/test_data_manager.cpp examples/data_manager_demo.cpp
mv examples/test_plasticity_verification.cpp examples/plasticity_verification.cpp
mv examples/test_plasticity_pcg.cpp examples/plasticity_pcg_demo.cpp
mv examples/test_plasticity_nr_complete.cpp examples/plasticity_nr_demo.cpp
mv examples/test_pure_shear.cpp examples/pure_shear_demo.cpp
mv examples/test_biaxial_tension.cpp examples/biaxial_tension_demo.cpp
mv examples/test_incremental_postprocess.cpp examples/incremental_postprocess_demo.cpp
mv examples/test_reaction_forces.cpp examples/reaction_forces_demo.cpp
```

### Step 3: 更新 CMakeLists.txt
需要移除已删除文件的编译目标，更新重命名文件的目标名称。

### Step 4: 运行测试确保功能完整
```bash
cd build
cmake ..
make
./bin/fem_tests
```

## 预期效果

- 删除文件：7 个
- 重命名文件：9 个  
- 保留核心示例：11 个
- 总计示例：20 个（从 27 个减少）
- 减少代码：~25%

## 验证清单

- [ ] 所有核心示例编译通过
- [ ] 所有测试通过
- [ ] 文档更新（README, TODO）
- [ ] Git 提交
