/**
 * poisson_2d_v2.cpp
 * 
 * 使用新 Assembler 接口的 2D Poisson 方程求解示例
 * 
 * 问题设定:
 * - 区域: [0,1]x[0,1]
 * - u = 0 (Dirichlet 边界条件)
 * - f = 1 (源项)
 * - k = 1 (导热系数)
 * 
 * 流程: 网格生成 → Assembler装配 → 边界条件 → 求解
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "math/cg.h"
#include "core/timer.h"
#include "core/logger.h"

#include <cmath>

using namespace fem;

// 计算三角形面积
Real triangle_area(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
    Real x1 = p1[0] - p0[0], y1 = p1[1] - p0[1];
    Real x2 = p2[0] - p0[0], y2 = p2[1] - p0[1];
    return 0.5 * std::abs(x1*y2 - x2*y1);
}

// 计算形函数梯度 (2D 三角形)
void triangle_gradients(const Vec3& p0, const Vec3& p1, const Vec3& p2, 
                       Real grad[][2]) {
    Real area = triangle_area(p0, p1, p2);
    if (area < 1e-15) return;
    
    // ∇N0
    grad[0][0] = (p1[1] - p2[1]) / (2.0 * area);
    grad[0][1] = (p2[0] - p1[0]) / (2.0 * area);
    
    // ∇N1
    grad[1][0] = (p2[1] - p0[1]) / (2.0 * area);
    grad[1][1] = (p0[0] - p2[0]) / (2.0 * area);
    
    // ∇N2
    grad[2][0] = (p0[1] - p1[1]) / (2.0 * area);
    grad[2][1] = (p1[0] - p0[0]) / (2.0 * area);
}

int main() {
    FEM_INFO("=== 2D Poisson Equation Demo (v2 - new Assembler) ===");
    
    // ── 1. 创建模型 ──
    Model model("Poisson 2D");
    
    // 材料 (导热系数)
    int mat_id = model.add_material("Conductivity");
    model.material(mat_id).set_property("k", 1.0);
    
    // 网格 (2D 三角形网格)
    int mesh_id = model.add_mesh("domain", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 20x20 网格
    MeshGenerator::generate_unit_square_tri(20, 20, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // ── 2. 创建 Assembler ──
    Timer timer;
    timer.start();
    
    Assembler assembler(model, 1);  // 标量场 (1 DOF per node)
    
    // 单元矩阵/向量计算函数
    auto elem_func = [](Index elem_id, const Mesh& mesh, DenseMatrix& Ke, Vector& Fe) {
        const Element& elem = mesh.element(elem_id);
        
        if (elem.type() != ElementType::Tri3) return;
        
        // 获取节点坐标
        const auto& nodes = elem.nodes();
        Vec3 coords[3];
        for (int i = 0; i < 3; ++i) {
            coords[i] = mesh.node(nodes[i]).coords();
        }
        
        // 计算单元刚度矩阵 Ke[i][j] = ∫_Ω k ∇Ni · ∇Nj dΩ
        Real area = triangle_area(coords[0], coords[1], coords[2]);
        Real grad[3][2];
        triangle_gradients(coords[0], coords[1], coords[2], grad);
        
        Real k = mesh.material()->property("k", 1.0);
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Ke(i, j) = k * (grad[i][0]*grad[j][0] + grad[i][1]*grad[j][1]) * area;
            }
        }
        
        // 计算单元载荷向量 Fe[i] = ∫_Ω f * Ni dΩ (f=1)
        Real source = 1.0;
        for (int i = 0; i < 3; ++i) {
            Fe[i] = source * area / 3.0;
        }
    };
    
    // 装配
    assembler.assemble(elem_func);
    FEM_INFO("Assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 3. 施加边界条件 ──
    timer.start();
    
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 0.0});
    bcs.push_back({"right", 0, 0.0});
    bcs.push_back({"bottom", 0, 0.0});
    bcs.push_back({"top", 0, 0.0});
    
    assembler.apply_dirichlet(bcs);
    FEM_INFO("Boundary conditions applied in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 4. 获取系统矩阵和向量 ──
    timer.start();
    SparseMatrixCSR K = assembler.matrix();
    const Vector& F = assembler.rhs();
    
    std::vector<Real> F_std = F;
    std::vector<Real> u_std(F.size(), 0.0);
    
    FEM_INFO("System retrieved in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 5. 求解线性系统 Ku = F ──
    timer.start();
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(1000);
    
    auto result = solver.solve(K, F_std, u_std);
    
    FEM_INFO("Solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("CG convergence: " + std::string(result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(result.residual));
    
    if (result.converged) {
        // ── 6. 输出结果 ──
        Real max_u = 0.0;
        std::size_t max_idx = 0;
        for (std::size_t i = 0; i < u_std.size(); ++i) {
            if (std::abs(u_std[i]) > max_u) {
                max_u = std::abs(u_std[i]);
                max_idx = i;
            }
        }
        
        const auto& max_coord = mesh.node(max_idx).coords();
        FEM_INFO("Max solution: u = " + fmt_sci(max_u) + 
                 " at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
        
        FEM_INFO("Demo completed successfully!");
    } else {
        FEM_ERROR("Solver failed to converge");
        return 1;
    }
    
    return 0;
}
