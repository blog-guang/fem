# ══════════════════════════════════════════════════════════════════
# FEM 示例程序
# ══════════════════════════════════════════════════════════════════

# ── 通用宏：减少重复代码 ──
macro(add_fem_example name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE fem)
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
endmacro()

# ══════════════════════════════════════════════════════════════════
# 基础示例（单物理场）
# ══════════════════════════════════════════════════════════════════

# Poisson 方程
add_fem_example(poisson_2d_v2)

# 热传导
add_fem_example(heat_2d)

# 线弹性
add_fem_example(elasticity_2d)

# 悬臂梁（Neumann BC 演示）- 待更新到新 API
# add_fem_example(cantilever_beam)

# ══════════════════════════════════════════════════════════════════
# 多物理场耦合
# ══════════════════════════════════════════════════════════════════

# 热-结构耦合 - 待更新到新 API
# add_fem_example(thermal_stress_2d)

# ══════════════════════════════════════════════════════════════════
# 非线性问题
# ══════════════════════════════════════════════════════════════════

# 几何非线性桁架
add_fem_example(nonlinear_truss)

# 材料非线性（塑性）
add_fem_example(plasticity_demo)
add_fem_example(plasticity_pcg_demo)
add_fem_example(plasticity_nr_demo)
add_fem_example(plasticity_verification)

# 多轴应力状态
add_fem_example(pure_shear_demo)
add_fem_example(biaxial_tension_demo)

# ══════════════════════════════════════════════════════════════════
# 材料模型演示
# ══════════════════════════════════════════════════════════════════

add_fem_example(material_demo)
add_fem_example(material_simple_test)

# ══════════════════════════════════════════════════════════════════
# 数据管理和后处理
# ══════════════════════════════════════════════════════════════════

add_fem_example(data_manager_demo)
add_fem_example(incremental_postprocess_demo)
add_fem_example(reaction_forces_demo)

# ══════════════════════════════════════════════════════════════════
# 性能测试和对比
# ══════════════════════════════════════════════════════════════════

add_fem_example(benchmark_preconditioners)
