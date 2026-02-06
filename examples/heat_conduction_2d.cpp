/**
 * heat_conduction_2d.cpp
 *
 * 二维稳态热传导方程求解:
 *   -∇·(k∇T) = Q,   Ω = [0,1] x [0,1]
 *         T   = 0,   ∂Ω (所有边界 Dirichlet, 固定温度)
 *
 * 流程: 生成网格 → 装配 K+F → 施加 BC → CG 求解 → VTK 输出
 */

#include "core/logger.h"
#include "core/timer.h"
#include "mesh/mesh_generator.h"
#include "element/element_base.h"
#include "assembly/assembler.h"
#include "assembly/sparse_matrix.h"
#include "assembly/boundary_condition.h"
#include "solver/solver.h"
#include "io/io.h"
#include "physics/heat_conduction.h"

#include <cmath>
#include <iostream>

using namespace fem;

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    Timer total;
    total.start();

    // ── 1. 生成网格 ──
    constexpr std::size_t N = 20;  // NxN 划分
    Mesh mesh = generate_unit_square_tri(N, N);

    // ── 2. 设置材料参数 ──
    HeatMaterial material;
    material.conductivity = 1.0;  // k = 1
    material.source = 1.0;        // Q = 1

    // ── 3. 装配 ──
    auto elem = create_element(ElementType::Triangle2D);
    Assembler assembler(mesh, *elem, 1);  // 1 DOF per node for scalar field

    COOMatrix K_coo;
    Vector    F;
    assembler.assemble(heat_stiffness, heat_load, K_coo, F, &material);
    FEM_INFO("Assembly done. nnz(COO) = " + std::to_string(K_coo.values.size()));

    // COO → CSR
    CSRMatrix K = coo_to_csr(K_coo);
    FEM_INFO("CSR nnz = " + std::to_string(K.values.size()));

    // ── 4. 边界条件: 四边 T=0 ──
    std::vector<BoundaryCondition> bcs = {
        {BCType::Dirichlet, "bottom", 0.0, -1},  // -1 means all components (scalar)
        {BCType::Dirichlet, "top",    0.0, -1},
        {BCType::Dirichlet, "left",   0.0, -1},
        {BCType::Dirichlet, "right",  0.0, -1}
    };
    apply_boundary_conditions(K, F, mesh, bcs, 1);  // 1 DOF per node
    FEM_INFO("Boundary conditions applied.");

    // ── 5. 求解 ──
    auto solver = create_solver(SolverType::CG);
    solver->set_tol(1e-10);
    solver->set_max_iter(2000);

    Vector T;
    SolveResult result = solver->solve(K, F, T);

    FEM_INFO("Solve done: converged=" + std::to_string(result.converged) +
             " iters=" + std::to_string(result.iterations) +
             " residual=" + fmt_sci(result.residual));

    // ── 6. 输出 ──
    VTKWriter writer("./output");
    writer.write("heat_conduction_2d.vtk", mesh, {{"temperature", &T}});

    // ── 7. 简单验证: 中心点处 T 的值 ──
    // 找最接近 (0.5, 0.5) 的节点
    Real best_dist = 1e30;
    Index best_idx = 0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real dx = mesh.coords(i)[0] - 0.5;
        Real dy = mesh.coords(i)[1] - 0.5;
        Real d  = dx*dx + dy*dy;
        if (d < best_dist) { best_dist = d; best_idx = i; }
    }
    FEM_INFO("T at center ≈ " + std::to_string(T[best_idx]));

    total.stop();
    Timer::print("Total", total.elapsed_s());

    return 0;
}