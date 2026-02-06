/**
 * thermal_stress_2d.cpp
 *
 * 二维热-结构耦合问题 (单向耦合):
 *   步骤1: 求解热传导 → 温度场 T
 *   步骤2: 用温度场计算热应变载荷 → 求解弹性
 *
 * 流程: 热传导求解 → 热应变计算 → 结构分析 → VTK 输出
 */

#include "core/logger.h"
#include "core/timer.h"
#include "mesh/mesh_generator.h"
#include "element/element_base.h"
#include "assembly/assembler.h"
#include <algorithm>
#include "assembly/sparse_matrix.h"
#include "assembly/boundary_condition.h"
#include "solver/solver.h"
#include "io/io.h"
#include "physics/heat_conduction.h"
#include "physics/elasticity.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace fem;

// 扩展ElasticCtx来包含温度信息
struct ThermalElasticCtx : public ElasticCtx {
    std::vector<Real> temperature;  // 单元温度
    Real alpha{1e-5};              // 热膨胀系数
};

// 扩展的弹性载荷函数，包括热应变
void thermal_elasticity_load(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* ctx) {
    // 先计算体力载荷
    elasticity_load(coords, n_nodes, dofs_per_node, Fe, ctx);
    
    // 添加热应变载荷
    const ThermalElasticCtx* thermal_ctx = static_cast<ThermalElasticCtx*>(ctx);
    const ElasticMaterial& mat = thermal_ctx->mat;
    Real alpha = thermal_ctx->alpha;
    
    if (n_nodes != 3 || dofs_per_node != 2) {
        return;
    }

    Real E = mat.E;
    Real nu = mat.poisson;
    Real t = mat.thickness;
    
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    Real area = std::abs(detJ) * 0.5;

    // 计算平均温度作为单元温度
    Real avg_temp = 0.0;
    for (std::size_t i = 0; i < n_nodes; ++i) {
        avg_temp += thermal_ctx->temperature[i];
    }
    avg_temp /= n_nodes;
    
    // 计算热应变: ε_thermal = [alpha*T, alpha*T, 0]
    Real eps_xx_th = alpha * avg_temp;
    Real eps_yy_th = alpha * avg_temp;
    
    // B矩阵 (重新计算)
    Real b1 = y1 - y2;
    Real b2 = y2 - y0;
    Real b3 = y0 - y1;
    Real c1 = x2 - x1;
    Real c2 = x0 - x2;
    Real c3 = x1 - x0;

    Real factor = 1.0 / (2.0 * area);
    
    Real B[3][6] = {
        {factor * b1, 0.0, factor * b2, 0.0, factor * b3, 0.0},      // strain_x = du/dx
        {0.0, factor * c1, 0.0, factor * c2, 0.0, factor * c3},      // strain_y = dv/dy
        {factor * c1, factor * b1, factor * c2, factor * b2, factor * c3, factor * b3}  // gamma_xy = du/dy + dv/dx
    };

    // D矩阵 (平面应力)
    Real D[3][3];
    Real D_factor = E / (1.0 - nu * nu);
    D[0][0] = D_factor;           D[0][1] = D_factor * nu;    D[0][2] = 0.0;
    D[1][0] = D_factor * nu;      D[1][1] = D_factor;         D[1][1] = D_factor;
    D[2][0] = 0.0;                D[2][1] = 0.0;              D[2][2] = D_factor * (1.0 - nu) / 2.0;

    // 计算热应变向量
    Real eps_th[3] = {eps_xx_th, eps_yy_th, 0.0};
    Real sigma_th[3];  // 热应力
    
    // sigma_th = D * eps_th
    for (int i = 0; i < 3; ++i) {
        sigma_th[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            sigma_th[i] += D[i][j] * eps_th[j];
        }
    }
    
    // 计算热应变等效节点力: Fe_thermal = t * area * B^T * sigma_th
    Real Fe_thermal[6];
    for (int i = 0; i < 6; ++i) {
        Fe_thermal[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            Fe_thermal[i] += B[j][i] * sigma_th[j];
        }
        Fe_thermal[i] *= t * area;
    }

    // 将热应变载荷添加到总载荷
    for (int i = 0; i < 6; ++i) {
        Fe[i] += Fe_thermal[i];
    }
}

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    Timer total;
    total.start();

    // ── 1. 生成网格 ──
    constexpr std::size_t N = 20;  // NxN 划分
    Mesh mesh = generate_unit_square_tri(N, N);

    // ── 2. 步骤1: 求解热传导问题 ──
    FEM_INFO("Step 1: Solving heat conduction...");
    
    // 设置热传导材料参数
    HeatMaterial heat_mat;
    heat_mat.conductivity = 1.0;
    heat_mat.source = 1.0;

    // 装配热传导系统
    auto elem = create_element(ElementType::Triangle2D);
    Assembler heat_assembler(mesh, *elem, 1);  // 1 DOF per node for temperature

    COOMatrix K_heat_coo;
    Vector    F_heat;
    heat_assembler.assemble(heat_stiffness, heat_load, K_heat_coo, F_heat, &heat_mat);
    
    // COO → CSR
    CSRMatrix K_heat = coo_to_csr(K_heat_coo);

    // 热传导边界条件: 四边温度为0
    std::vector<BoundaryCondition> heat_bcs = {
        {BCType::Dirichlet, "bottom", 0.0, -1},
        {BCType::Dirichlet, "top",    0.0, -1},
        {BCType::Dirichlet, "left",   0.0, -1},
        {BCType::Dirichlet, "right",  0.0, -1}
    };
    apply_boundary_conditions(K_heat, F_heat, mesh, heat_bcs, 1);

    // 求解温度场
    auto solver = create_solver(SolverType::CG);
    solver->set_tol(1e-10);
    solver->set_max_iter(2000);

    Vector T;
    SolveResult heat_result = solver->solve(K_heat, F_heat, T);
    FEM_INFO("Heat solve: converged=" + std::to_string(heat_result.converged) +
             " iters=" + std::to_string(heat_result.iterations));

    // ── 3. 步骤2: 求解热弹性问题 ──
    FEM_INFO("Step 2: Solving thermal stress...");
    
    // 设置弹性材料参数
    ThermalElasticCtx elastic_ctx;
    elastic_ctx.mat.E = 1e6;
    elastic_ctx.mat.poisson = 0.3;
    elastic_ctx.mat.thickness = 1.0;
    elastic_ctx.load.fx = 0.0;
    elastic_ctx.load.fy = 0.0;  // 无体力
    elastic_ctx.alpha = 1e-5;   // 热膨胀系数
    
    // 为每个单元计算平均温度
    elastic_ctx.temperature.resize(mesh.num_cells());
    for (std::size_t c = 0; c < mesh.num_cells(); ++c) {
        const Cell& cell = mesh.cell(c);
        Real avg_temp = 0.0;
        for (std::size_t i = 0; i < cell.num_nodes; ++i) {
            avg_temp += T[cell.node(i)];
        }
        avg_temp /= cell.num_nodes;
        elastic_ctx.temperature[c] = avg_temp;
    }

    // 装配弹性系统（带热应变）
    Assembler elastic_assembler(mesh, *elem, 2);  // 2 DOFs per node (ux, uy)

    COOMatrix K_elastic_coo;
    Vector    F_elastic;
    elastic_assembler.assemble(elasticity_stiffness, thermal_elasticity_load, K_elastic_coo, F_elastic, &elastic_ctx);
    
    // COO → CSR
    CSRMatrix K_elastic = coo_to_csr(K_elastic_coo);

    // 弹性边界条件: 左边固定 ux=uy=0
    std::vector<BoundaryCondition> elastic_bcs = {
        {BCType::Dirichlet, "left", 0.0, 0},  // x方向位移为0
        {BCType::Dirichlet, "left", 0.0, 1}   // y方向位移为0
    };
    apply_boundary_conditions(K_elastic, F_elastic, mesh, elastic_bcs, 2);

    // 求解位移场
    Vector u;
    SolveResult elastic_result = solver->solve(K_elastic, F_elastic, u);
    FEM_INFO("Elastic solve: converged=" + std::to_string(elastic_result.converged) +
             " iters=" + std::to_string(elastic_result.iterations));

    // ── 4. 输出 ──
    // 温度场
    VTKWriter writer("./output");
    writer.write("thermal_stress_temperature.vtk", mesh, {{"temperature", &T}});
    
    // 位移场
    Vector ux(mesh.num_nodes()), uy(mesh.num_nodes());
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        ux[i] = u[2*i];     // x方向位移
        uy[i] = u[2*i+1];   // y方向位移
    }
    writer.write("thermal_stress_displacement.vtk", mesh, {{"ux", &ux}, {"uy", &uy}});

    // ── 5. 简单验证 ──
    FEM_INFO("Max temperature: " + std::to_string(*std::max_element(T.begin(), T.end())));
    FEM_INFO("Max displacement: ux=" + std::to_string(*std::max_element(ux.begin(), ux.end())) + 
             ", uy=" + std::to_string(*std::min_element(uy.begin(), uy.end())));

    total.stop();
    Timer::print("Total", total.elapsed_s());

    return 0;
}