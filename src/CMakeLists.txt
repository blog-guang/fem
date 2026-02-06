# ── fem 静态库 ──
# 按模块逐步添加，当前: Core

add_library(fem STATIC
    # Core
    core/logger.cpp
    core/timer.cpp
    # Mesh
    mesh/mesh.cpp
    mesh/mesh_generator.cpp
    # Element
    element/element_base.cpp
    element/triangle2d.cpp
    element/quad2d.cpp
    # Assembly
    assembly/sparse_matrix.cpp
    assembly/assembler.cpp
    assembly/boundary_condition.cpp
    # Solver
    solver/solver.cpp
    solver/cg.cpp
    solver/bicgstab.cpp
    # IO
    io/io.cpp
    # Physics
    physics/heat_conduction.cpp
    physics/elasticity.cpp
)

target_include_directories(fem PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

target_compile_features(fem PUBLIC cxx_std_17)
