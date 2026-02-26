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
