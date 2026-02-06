/**
 * elasticity_simple.cpp - 最小弹性力学测试
 * 
 * 2x2 网格，完全约束边界条件
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_v2.h"
#include "solver/cg.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== Minimal Elasticity Test ===");
    
    // 创建模型
    Model model("test");
    int mat_id = model.add_material("mat");
    model.material(mat_id).set_property("E", 1000.0);
    model.material(mat_id).set_property("nu", 0.3);
    
    int mesh_id = model.add_mesh("mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 最小网格: 2x2
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // 创建物理模块
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    // 装配
    Assembler assembler(model, 2);
    
    auto elem_func = [&elast](Index elem_id, const Mesh& mesh,
                              DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    assembler.assemble(elem_func);
    FEM_INFO("Assembly: " + std::to_string(assembler.num_dofs()) + " DOFs");
    
    // 边界条件：所有边界固定
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 0.0});
    bcs.push_back({"left", 1, 0.0});
    bcs.push_back({"right", 0, 0.0});
    bcs.push_back({"right", 1, 0.0});
    bcs.push_back({"bottom", 0, 0.0});
    bcs.push_back({"bottom", 1, 0.0});
    bcs.push_back({"top", 0, 0.0});
    bcs.push_back({"top", 1, 0.0});
    
    assembler.apply_dirichlet(bcs);
    
    // 获取系统
    SparseMatrixCSR K = assembler.matrix();
    const Vector& F = assembler.rhs();
    
    FEM_INFO("Matrix: " + std::to_string(K.nnz()) + " non-zeros");
    
    // 检查矩阵
    bool has_nan = false;
    for (std::size_t i = 0; i < K.values().size(); ++i) {
        if (std::isnan(K.values()[i]) || std::isinf(K.values()[i])) {
            FEM_ERROR("Matrix has NaN/Inf at index " + std::to_string(i));
            has_nan = true;
            break;
        }
    }
    
    if (has_nan) {
        return 1;
    }
    
    // 求解
    std::vector<Real> F_std = F.raw();
    std::vector<Real> u_std(F.size(), 0.0);
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(100);
    
    auto result = solver.solve(K, F_std, u_std);
    
    FEM_INFO("CG: " + std::string(result.converged ? "converged" : "failed") + 
             ", residual=" + fmt_sci(result.residual));
    
    if (result.converged) {
        FEM_INFO("All boundaries fixed → solution should be near zero");
        Real max_u = 0.0;
        for (std::size_t i = 0; i < u_std.size(); ++i) {
            max_u = std::max(max_u, std::abs(u_std[i]));
        }
        FEM_INFO("Max |u| = " + fmt_sci(max_u));
        
        if (max_u < 1e-6) {
            FEM_INFO("Test PASSED!");
            return 0;
        } else {
            FEM_ERROR("Test FAILED: expected near-zero solution");
            return 1;
        }
    } else {
        FEM_ERROR("Solver failed");
        return 1;
    }
}
