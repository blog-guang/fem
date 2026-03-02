# ── GoogleTest (submodule) ──
add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/googletest 
                 ${CMAKE_BINARY_DIR}/third_party/googletest
                 EXCLUDE_FROM_ALL)

# ── 测试目标 (按模块逐步添加) ──
add_executable(fem_tests
    test_core.cpp
    test_math.cpp
    test_material.cpp
    test_material_factory.cpp
    test_material_models.cpp
    test_advanced_materials.cpp
    test_shape_function.cpp
    test_shape_integration.cpp
    test_solver.cpp
    test_solver_factory.cpp
    test_assembler.cpp
    test_physics.cpp
    test_physics_factory.cpp
    test_mesh_v2.cpp
    test_mesh_generator.cpp
    test_io.cpp
    test_newton_raphson.cpp
    test_data.cpp
    test_ilu.cpp
    test_amg.cpp
    test_load_step.cpp
    test_neo_hookean.cpp
)

target_link_libraries(fem_tests PRIVATE fem gtest_main)
target_include_directories(fem_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME fem_tests COMMAND fem_tests)
