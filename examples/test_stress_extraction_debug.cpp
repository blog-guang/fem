/**
 * test_stress_extraction_debug.cpp - 调试应力提取逻辑
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "solver/pcg.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

int main() {
    Logger::instance().set_level(LogLevel::WARN);
    
    std::cout << "=== Stress Extraction Debug ===\n\n";
    
    // 几何参数（一致单位：mm）
    Real length = 10.0;   // mm
    Real width = 1.0;     // mm
    
    // 材料参数
    Real E = 200e3;       // MPa = N/mm²
    Real nu = 0.3;
    
    std::cout << "Geometry:\n";
    std::cout << "  Length: " << length << " mm\n";
    std::cout << "  Width:  " << width << " mm\n";
    std::cout << "  Area:   " << width << " mm²\n\n";
    
    std::cout << "Material:\n";
    std::cout << "  E = " << E << " MPa\n";
    std::cout << "  ν = " << nu << "\n\n";
    
    // 创建网格
    Model model("debug");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("specimen", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    int nx = 4, ny = 2;
    MeshGenerator::generate_unit_square_quad(nx, ny, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    // 缩放
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Vec3& coords = mesh.node(i).coords();
        coords[0] *= length;
        coords[1] *= width;
    }
    
    std::cout << "Mesh:\n";
    std::cout << "  Nodes:    " << mesh.num_nodes() << "\n";
    std::cout << "  Elements: " << mesh.num_elements() << "\n\n";
    
    // 打印节点坐标
    std::cout << "Node Coordinates:\n";
    std::cout << "ID     X (mm)    Y (mm)\n";
    std::cout << "-------------------------\n";
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        const Vec3& coords = mesh.node(i).coords();
        std::cout << std::setw(2) << i << "  "
                  << std::fixed << std::setprecision(2)
                  << std::setw(8) << coords[0] << "  "
                  << std::setw(8) << coords[1] << "\n";
    }
    std::cout << "\n";
    
    // 识别边界节点
    std::vector<Index> left_nodes, right_nodes;
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Real x = mesh.node(i).coords()[0];
        if (std::abs(x) < 1e-6) {
            left_nodes.push_back(i);
        }
        if (std::abs(x - length) < 1e-6) {
            right_nodes.push_back(i);
        }
    }
    
    std::cout << "Boundary Nodes:\n";
    std::cout << "  Left:  " << left_nodes.size() << " nodes: ";
    for (auto id : left_nodes) std::cout << id << " ";
    std::cout << "\n  Right: " << right_nodes.size() << " nodes: ";
    for (auto id : right_nodes) std::cout << id << " ";
    std::cout << "\n\n";
    
    // 创建材料和物理模块
    IsotropicElastic material(E, nu, 2, true);  // 2D plane_stress
    ElasticityUnified physics(&material, 2);
    
    // 装配
    Assembler assembler(model, 2);
    assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        physics.compute_element(elem_id, m, Ke, Fe);
    });
    
    // 边界条件
    Real displacement = 0.001;  // mm (0.1% 应变)
    
    std::vector<DirichletBC> bcs;
    for (auto id : left_nodes) {
        bcs.push_back({"left", static_cast<Index>(id * 2), 0.0});      // u_x = 0
    }
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        if (std::abs(mesh.node(i).coords()[1]) < 1e-6) {
            bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});  // u_y = 0
        }
    }
    for (auto id : right_nodes) {
        bcs.push_back({"right", static_cast<Index>(id * 2), displacement});  // u_x = u
    }
    
    assembler.apply_dirichlet(bcs);
    
    std::cout << "Applied " << bcs.size() << " boundary conditions\n";
    std::cout << "  Displacement: " << displacement << " mm\n";
    std::cout << "  Strain:       " << displacement/length << "\n\n";
    
    // 求解
    PCGSolver solver("jacobi");
    solver.set_max_iter(10000);
    solver.set_tol(1e-10);
    
    const auto& K = assembler.matrix();
    const auto& F = assembler.rhs();
    std::vector<Real> u(F.size(), 0.0);
    
    auto result = solver.solve(K, F.raw(), u);
    
    if (!result.converged) {
        std::cerr << "Solver failed!\n";
        return 1;
    }
    
    std::cout << "Solver converged in " << result.iterations << " iterations\n\n";
    
    // 打印位移
    std::cout << "Displacement Solution:\n";
    std::cout << "Node   u_x (mm)   u_y (mm)\n";
    std::cout << "-----------------------------\n";
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        std::cout << std::setw(4) << i << "  "
                  << std::fixed << std::setprecision(6)
                  << std::setw(10) << u[i*2] << "  "
                  << std::setw(10) << u[i*2+1] << "\n";
    }
    std::cout << "\n";
    
    // 方法1：从内力计算应力（通过单元积分）
    // 对于单轴拉伸，应力应该是均匀的
    // σ = E * ε = E * (u/L)
    
    Real avg_displacement = 0.0;
    for (auto id : right_nodes) {
        avg_displacement += u[id * 2];
    }
    avg_displacement /= right_nodes.size();
    
    Real strain_from_disp = avg_displacement / length;
    Real stress_from_strain = E * strain_from_disp;
    
    std::cout << "Method 1: From Displacement\n";
    std::cout << "  Average displacement (right): " << avg_displacement << " mm\n";
    std::cout << "  Strain:  ε = " << strain_from_disp << "\n";
    std::cout << "  Stress:  σ = E*ε = " << stress_from_strain << " MPa\n\n";
    
    // 方法2：从反力计算（需要正确处理）
    // 反力 = K*u - F_external
    // 但在位移边界条件下，F_external = 0
    // 所以反力 = K*u（在约束节点上）
    
    std::vector<Real> Ku(u.size(), 0.0);
    K.matvec(u.data(), Ku.data());
    
    // 从 F_external 计算反力
    // reaction[i] = Ku[i] - F[i]
    std::vector<Real> reaction(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
        reaction[i] = Ku[i] - F.raw()[i];
    }
    
    std::cout << "Method 2: From Reaction Forces\n";
    std::cout << "Left boundary reactions (x-direction):\n";
    Real F_left_total = 0.0;
    for (auto id : left_nodes) {
        Real fx = reaction[id * 2];
        F_left_total += fx;
        std::cout << "  Node " << id << ": " << fx << " N\n";
    }
    std::cout << "  Total: " << F_left_total << " N\n\n";
    
    // 计算应力
    Real area = width * 1.0;  // mm² (假设单位厚度)
    Real stress_from_reaction = std::abs(F_left_total) / area;
    
    std::cout << "=== Stress Calculation ===\n";
    std::cout << "Cross-sectional area: " << area << " mm²\n";
    std::cout << "From displacement:  σ = " << stress_from_strain << " MPa\n";
    std::cout << "From reaction:      σ = |" << F_left_total << "| / " << area << " = " << stress_from_reaction << " MPa\n\n";
    
    // 解析解
    Real strain = displacement / length;
    Real sigma_theory = E * strain;
    
    std::cout << "=== Analytical Solution ===\n";
    std::cout << "Strain:  ε = " << strain << "\n";
    std::cout << "Stress:  σ = E*ε = " << E << " * " << strain << " = " << sigma_theory << " MPa\n\n";
    
    // 误差
    Real error_disp = std::abs(stress_from_strain - sigma_theory) / sigma_theory * 100.0;
    Real error_reaction = std::abs(stress_from_reaction - sigma_theory) / sigma_theory * 100.0;
    
    std::cout << "=== Error Analysis ===\n";
    std::cout << "From displacement: " << error_disp << "%\n";
    std::cout << "From reaction:     " << error_reaction << "%\n\n";
    
    if (error_disp < 5.0) {
        std::cout << "✓ Displacement-based stress extraction is accurate!\n";
    } else {
        std::cout << "✗ Displacement-based stress has error: " << error_disp << "%\n";
    }
    
    if (error_reaction < 5.0) {
        std::cout << "✓ Reaction-based stress extraction is accurate!\n";
    } else {
        std::cout << "✗ Reaction-based stress has error: " << error_reaction << "%\n";
        std::cout << "   (This may be expected due to boundary condition handling)\n";
    }
    
    return 0;
}
