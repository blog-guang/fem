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
    # Material (本构模型)
    material/state_variables.cpp
    material/material.cpp
    material/isotropic_elastic.cpp
    material/orthotropic_elastic.cpp
    material/j2_plasticity.cpp
    material/j2_plasticity_kinematic.cpp
    # Shape (形函数)
    shape/shape_function.cpp
    shape/shape_tri3.cpp
    shape/shape_quad4.cpp
    shape/shape_tet4.cpp
    shape/shape_brick8.cpp
    shape/shape_function_factory.cpp
    # Mesh (新架构)
    mesh/material.cpp
    mesh/mesh.cpp
    mesh/model.cpp
    mesh/mesh_generator.cpp
    # Assembly
    assembly/assembler.cpp
    # Physics
    physics/physics_base.cpp
    physics/heat_unified.cpp
    physics/elasticity_unified.cpp
    # Data (场数据管理)
    data/base_data.cpp
    # PostProcess (后处理)
    postprocess/post_processor.cpp
    postprocess/post_processor_incremental.cpp
    # IO
    io/vtk_writer.cpp
    # Math - Linear Solvers
    math/solver.cpp
    math/cg.cpp
    math/pcg.cpp
    math/bicgstab.cpp
    math/newton_raphson.cpp
)

target_include_directories(fem PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

target_compile_features(fem PUBLIC cxx_std_17)
