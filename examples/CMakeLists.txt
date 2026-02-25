# ── 示例程序 ──

# 材料模型演示
add_executable(material_demo material_demo.cpp)
target_link_libraries(material_demo PRIVATE fem)
target_include_directories(material_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 材料模型演示
add_executable(material_demo material_demo.cpp)
target_link_libraries(material_demo PRIVATE fem)
target_include_directories(material_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 材料模型简化测试
add_executable(material_simple_test material_simple_test.cpp)
target_link_libraries(material_simple_test PRIVATE fem)
target_include_directories(material_simple_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
