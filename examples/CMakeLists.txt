# ── 示例程序 ──

# 材料模型演示
add_executable(material_demo material_demo.cpp)
target_link_libraries(material_demo PRIVATE fem)
target_include_directories(material_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 材料模型简化测试
add_executable(material_simple_test material_simple_test.cpp)
target_link_libraries(material_simple_test PRIVATE fem)
target_include_directories(material_simple_test PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 悬臂梁材料本构分析（待修复）
# add_executable(cantilever_material cantilever_material.cpp)
# target_link_libraries(cantilever_material PRIVATE fem)
# target_include_directories(cantilever_material PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 统一物理模块测试
add_executable(test_unified_heat test_unified_heat.cpp)
target_link_libraries(test_unified_heat PRIVATE fem)
target_include_directories(test_unified_heat PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(test_unified_elasticity test_unified_elasticity.cpp)
target_link_libraries(test_unified_elasticity PRIVATE fem)
target_include_directories(test_unified_elasticity PRIVATE ${CMAKE_SOURCE_DIR}/src)

# ElasticityUnified 与自定义材料集成示例
add_executable(test_elasticity_with_plasticity test_elasticity_with_plasticity.cpp)
target_link_libraries(test_elasticity_with_plasticity PRIVATE fem)
target_include_directories(test_elasticity_with_plasticity PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 数据管理系统示例
add_executable(test_data_manager test_data_manager.cpp)
target_link_libraries(test_data_manager PRIVATE fem)
target_include_directories(test_data_manager PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 弹塑性测试（材料非线性）
# 应力提取调试
add_executable(test_stress_extraction_debug test_stress_extraction_debug.cpp)
target_link_libraries(test_stress_extraction_debug PRIVATE fem)
target_include_directories(test_stress_extraction_debug PRIVATE ${CMAKE_SOURCE_DIR}/src)

# PCG 求解器测试（完整 FEM + 塑性）
add_executable(test_plasticity_pcg test_plasticity_pcg.cpp)
target_link_libraries(test_plasticity_pcg PRIVATE fem)
target_include_directories(test_plasticity_pcg PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Newton-Raphson 非线性求解（完整）
add_executable(test_plasticity_nr_complete test_plasticity_nr_complete.cpp)
target_link_libraries(test_plasticity_nr_complete PRIVATE fem)
target_include_directories(test_plasticity_nr_complete PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 纯剪切测试
add_executable(test_pure_shear test_pure_shear.cpp)
target_link_libraries(test_pure_shear PRIVATE fem)
target_include_directories(test_pure_shear PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 双轴拉伸测试
add_executable(test_biaxial_tension test_biaxial_tension.cpp)
target_link_libraries(test_biaxial_tension PRIVATE fem)
target_include_directories(test_biaxial_tension PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 增量后处理测试
add_executable(test_incremental_postprocess test_incremental_postprocess.cpp)
target_link_libraries(test_incremental_postprocess PRIVATE fem)
target_include_directories(test_incremental_postprocess PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 验证测试（单元级别）
add_executable(test_plasticity_verification test_plasticity_verification.cpp)
target_link_libraries(test_plasticity_verification PRIVATE fem)
target_include_directories(test_plasticity_verification PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 简单拉伸测试（待修复）
# add_executable(test_plasticity_simple test_plasticity_simple.cpp)
# target_link_libraries(test_plasticity_simple PRIVATE fem)
# target_include_directories(test_plasticity_simple PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 弯曲梁测试（待完善）
# add_executable(test_plasticity_bending test_plasticity_bending.cpp)
# target_link_libraries(test_plasticity_bending PRIVATE fem)
# target_include_directories(test_plasticity_bending PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 反力提取测试
add_executable(test_reaction_forces test_reaction_forces.cpp)
target_link_libraries(test_reaction_forces PRIVATE fem)
target_include_directories(test_reaction_forces PRIVATE ${CMAKE_SOURCE_DIR}/src)
