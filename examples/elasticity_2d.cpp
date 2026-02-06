/**
 * elasticity_2d.cpp
 *
 * 二维平面应力弹性问题求解 (悬臂梁):
 *   div(σ) + f = 0,   Ω = [0,2] x [0,1]
 *   σ = Dε, ε = sym(grad(u))
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
#include "physics/elasticity.h"

#include <cmath>
#include <iostream>

using namespace fem;

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    Timer total;
    total.start();

    // ── 1. 生成网格 ──
    constexpr std::size_t N = 40;  // NxN 划分，比例调整为 2:1 (矩形域 [0,2]x[0,1])
    Mesh mesh = generate_unit_square_tri(N, N/2);  // 生成 [0,1]x[0,0.5]，然后缩放到 [0,2]x[0,1]

    // 缩放网格到 [0,2]x[0,1]
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Vec3& coord = mesh.coords(i);
        coord[0] *= 2.0;  // x坐标乘以2
    }

    // ── 2. 设置材料参数 ──
    ElasticCtx ctx;
    ctx.mat.E = 1e6;          // 弹性模量
    ctx.mat.poisson = 0.3;    // 泊松比
    ctx.mat.thickness = 1.0;  // 厚度
    ctx.load.fx = 0.0;        // x方向体力
    ctx.load.fy = -1.0;       // y方向体力 (重力)

    // ── 3. 装配 ──
    auto elem = create_element(ElementType::Triangle2D);
    Assembler assembler(mesh, *elem, 2);  // 2 DOFs per node (ux, uy)

    COOMatrix K_coo;
    Vector    F;
    assembler.assemble(elasticity_stiffness, elasticity_load, K_coo, F, &ctx);
    FEM_INFO("Assembly done. nnz(COO) = " + std::to_string(K_coo.values.size()));

    // COO → CSR
    CSRMatrix K = coo_to_csr(K_coo);
    FEM_INFO("CSR nnz = " + std::to_string(K.values.size()));

    // ── 4. 边界条件: 左边 (x=0) 固定 ux=uy=0 ──
    std::vector<BoundaryCondition> bcs = {
        {BCType::Dirichlet, "left", 0.0, 0},  // x方向位移为0
        {BCType::Dirichlet, "left", 0.0, 1}   // y方向位移为0
    };
    apply_boundary_conditions(K, F, mesh, bcs, 2);  // 2 DOFs per node
    FEM_INFO("Boundary conditions applied.");

    // ── 5. 求解 ──
    auto solver = create_solver(SolverType::CG);
    solver->set_tol(1e-8);
    solver->set_max_iter(2000);

    Vector u;
    SolveResult result = solver->solve(K, F, u);

    FEM_INFO("Solve done: converged=" + std::to_string(result.converged) +
             " iters=" + std::to_string(result.iterations) +
             " residual=" + fmt_sci(result.residual));

    // ── 6. 输出 ──
    // 分离x和y方向的位移
    Vector ux(mesh.num_nodes()), uy(mesh.num_nodes());
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        ux[i] = u[2*i];     // x方向位移
        uy[i] = u[2*i+1];   // y方向位移
    }

    VTKWriter writer("./output");
    writer.write("elasticity_2d.vtk", mesh, {{"ux", &ux}, {"uy", &uy}});

    // ── 7. 简单验证: 右端最大变形 ──
    // 找最接近 (2.0, 0.5) 的节点（自由端中心）
    Real best_dist = 1e30;
    Index best_idx = 0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real dx = mesh.coords(i)[0] - 2.0;
        Real dy = mesh.coords(i)[1] - 0.5;
        Real d  = dx*dx + dy*dy;
        if (d < best_dist) { best_dist = d; best_idx = i; }
    }
    FEM_INFO("Max displacement at free end: ux=" + std::to_string(ux[best_idx]) + 
             ", uy=" + std::to_string(uy[best_idx]));

    total.stop();
    Timer::print("Total", total.elapsed_s());

    return 0;
}