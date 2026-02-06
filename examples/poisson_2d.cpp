/**
 * poisson_2d.cpp
 *
 * 二维 Poisson 方程求解:
 *   -∇²u = 1,   Ω = [0,1] x [0,1]
 *    u   = 0,   ∂Ω (所有边界 Dirichlet)
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

#include <cmath>
#include <iostream>

using namespace fem;

// ── 线性三角形单元刚度矩阵 ──
// -∇²u 的弱形式: K_ij = ∫ ∇N_i · ∇N_j dA
//
// 对线性三角形: Jacobian 为常数, 梯度为常数
//   J = [[x1-x0, x2-x0],
//        [y1-y0, y2-y0]]
//   detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
//   area = |detJ| / 2
//   dN/dx = J^{-T} * dN/dxi
//   K_ij  = area * (dNi_dx*dNj_dx + dNi_dy*dNj_dy)
//
static void poisson_stiffness(const Vec3* coords, std::size_t /*n_nodes*/, std::size_t /*dofs_per_node*/, Real* Ke, void* /*ctx*/) {
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    Real area = std::abs(detJ) * 0.5;

    // 参考梯度 dN/dxi:  N0=(-1,-1), N1=(1,0), N2=(0,1)
    // 物理梯度 dN/dx = J^{-T} * dN/dxi
    //   J^{-T} = (1/detJ) * [[ y2-y0, -(y1-y0)],
    //                         [-(x2-x0),  x1-x0]]
    Real dNdxi[3][2] = {{-1, -1}, {1, 0}, {0, 1}};
    Real dNdx[3][2];   // [node][x,y]

    for (int i = 0; i < 3; ++i) {
        dNdx[i][0] = ( (y2-y0)*dNdxi[i][0] - (y1-y0)*dNdxi[i][1] ) / detJ;
        dNdx[i][1] = (-(x2-x0)*dNdxi[i][0] + (x1-x0)*dNdxi[i][1] ) / detJ;
    }

    // Ke[i*3+j] = area * dot(dNdx[i], dNdx[j])
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke[i*3+j] = area * (dNdx[i][0]*dNdx[j][0] + dNdx[i][1]*dNdx[j][1]);
        }
    }
}

// ── 载荷向量 (f=1) ──
// F_i = ∫ N_i * f dA = f * area / 3  (线性三角形, 每节点均分)
//
static void poisson_load(const Vec3* coords, std::size_t /*n_nodes*/, std::size_t /*dofs_per_node*/, Real* Fe, void* /*ctx*/) {
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    Real area = std::abs(detJ) * 0.5;

    Fe[0] = Fe[1] = Fe[2] = area / 3.0;
}

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    Timer total;
    total.start();

    // ── 1. 生成网格 ──
    constexpr std::size_t N = 20;  // NxN 划分
    Mesh mesh = generate_unit_square_tri(N, N);

    // ── 2. 装配 ──
    auto elem = create_element(ElementType::Triangle2D);
    Assembler assembler(mesh, *elem, 1);  // 1 DOF per node for scalar field

    COOMatrix K_coo;
    Vector    F;
    assembler.assemble(poisson_stiffness, poisson_load, K_coo, F, nullptr);
    FEM_INFO("Assembly done. nnz(COO) = " + std::to_string(K_coo.values.size()));

    // COO → CSR
    CSRMatrix K = coo_to_csr(K_coo);
    FEM_INFO("CSR nnz = " + std::to_string(K.values.size()));

    // ── 3. 边界条件: 四边 u=0 ──
    std::vector<BoundaryCondition> bcs = {
        {BCType::Dirichlet, "bottom", 0.0, -1},  // -1 means all components (scalar)
        {BCType::Dirichlet, "top",    0.0, -1},
        {BCType::Dirichlet, "left",   0.0, -1},
        {BCType::Dirichlet, "right",  0.0, -1}
    };
    apply_boundary_conditions(K, F, mesh, bcs, 1);  // 1 DOF per node
    FEM_INFO("Boundary conditions applied.");

    // ── 4. 求解 ──
    auto solver = create_solver(SolverType::CG);
    solver->set_tol(1e-10);
    solver->set_max_iter(2000);

    Vector u;
    SolveResult result = solver->solve(K, F, u);

    FEM_INFO("Solve done: converged=" + std::to_string(result.converged) +
             " iters=" + std::to_string(result.iterations) +
             " residual=" + fmt_sci(result.residual));

    // ── 5. 输出 ──
    VTKWriter writer("./output");
    writer.write("poisson_2d.vtk", mesh, {{"u", &u}});

    // ── 6. 简单验证: 中心点处 u 的近似值 ──
    // 解析解 (正方形 Poisson, f=1, u=0 BC) 中心点 u ≈ 0.0737
    // 找最接近 (0.5, 0.5) 的节点
    Real best_dist = 1e30;
    Index best_idx = 0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real dx = mesh.coords(i)[0] - 0.5;
        Real dy = mesh.coords(i)[1] - 0.5;
        Real d  = dx*dx + dy*dy;
        if (d < best_dist) { best_dist = d; best_idx = i; }
    }
    FEM_INFO("u at center ≈ " + std::to_string(u[best_idx]) +
             " (analytical ≈ 0.0737)");

    total.stop();
    Timer::print("Total", total.elapsed_s());

    return 0;
}
