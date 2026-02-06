# ── fem 静态库 ──
# 按模块逐步添加，当前: Core

add_library(fem STATIC
    # Core
    core/logger.cpp
    core/timer.cpp
    # Mesh (新架构)
    mesh/material.cpp
    mesh/mesh.cpp
    mesh/model.cpp
    mesh/mesh_generator.cpp
    # Assembly (保留 sparse_matrix 供 solver 使用)
    assembly/sparse_matrix.cpp
    # Solver
    solver/solver.cpp
    solver/cg.cpp
    solver/bicgstab.cpp
)

target_include_directories(fem PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

target_compile_features(fem PUBLIC cxx_std_17)
