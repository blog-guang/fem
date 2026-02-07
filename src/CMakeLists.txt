# ── fem 静态库 ──
# 按模块逐步添加，当前: Core

add_library(fem STATIC
    # Core
    core/logger.cpp
    core/timer.cpp
    # Math
    math/vector.cpp
    math/dense_matrix.cpp
    math/sparse_matrix.cpp
    # Mesh (新架构)
    mesh/material.cpp
    mesh/mesh.cpp
    mesh/model.cpp
    mesh/mesh_generator.cpp
    # Assembly
    assembly/assembler.cpp
    # Physics
    physics/heat.cpp
    physics/elasticity_v2.cpp
    # IO
    io/vtk_writer.cpp
    # Solver
    solver/solver.cpp
    solver/cg.cpp
    solver/bicgstab.cpp
)

target_include_directories(fem PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

target_compile_features(fem PUBLIC cxx_std_17)
