/**
 * elasticity_full_debug.cpp - 完整系统调试
 * 
 * 检查装配 → 边界条件 → 求解的完整流程
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_v2.h"
#include "solver/cg.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

void print_matrix_info(const SparseMatrixCSR& K, const std::string& name) {
    FEM_INFO("\n=== " + name + " ===");
    FEM_INFO("Size: " + std::to_string(K.rows()) + "x" + std::to_string(K.cols()));
    FEM_INFO("Non-zeros: " + std::to_string(K.nnz()));
    
    // 检查 NaN/Inf
    bool has_nan = false;
    for (std::size_t i = 0; i < K.values().size(); ++i) {
        if (std::isnan(K.values()[i]) || std::isinf(K.values()[i])) {
            FEM_ERROR("Matrix has NaN/Inf at value index " + std::to_string(i));
            has_nan = true;
            break;
        }
    }
    
    if (!has_nan) {
        FEM_INFO("No NaN/Inf ✓");
    }
    
    // 检查对角元
    FEM_INFO("\nDiagonal elements:");
    for (std::size_t i = 0; i < std::min(K.rows(), (std::size_t)10); ++i) {
        Real diag = 0.0;
        bool found = false;
        
        for (Index k = K.row_ptr()[i]; k < K.row_ptr()[i + 1]; ++k) {
            if (K.col_indices()[k] == i) {
                diag = K.values()[k];
                found = true;
                break;
            }
        }
        
        if (found) {
            FEM_INFO("  K(" + std::to_string(i) + "," + std::to_string(i) + ") = " + fmt_sci(diag));
        } else {
            FEM_ERROR("  K(" + std::to_string(i) + "," + std::to_string(i) + ") = MISSING!");
        }
    }
    
    if (K.rows() > 10) {
        FEM_INFO("  ... (showing first 10)");
    }
}

int main() {
    FEM_INFO("=== Full System Debug ===");
    
    // 最小系统: 2x2 网格
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // 物理模块
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    // 装配
    Assembler assembler(model, 2);  // 矢量场
    
    auto elem_func = [&elast](Index elem_id, const Mesh& mesh,
                              DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    assembler.assemble(elem_func);
    
    FEM_INFO("\n=== After Assembly ===");
    FEM_INFO("Total DOFs: " + std::to_string(assembler.num_dofs()));
    
    // 获取装配后的矩阵（边界条件前）
    SparseMatrixCSR K_before = assembler.matrix();
    const Vector& F_before = assembler.rhs();
    
    print_matrix_info(K_before, "Matrix Before BC");
    
    FEM_INFO("\nRHS before BC (first 10):");
    for (std::size_t i = 0; i < std::min(F_before.size(), (std::size_t)10); ++i) {
        FEM_INFO("  F[" + std::to_string(i) + "] = " + fmt_sci(F_before[i]));
    }
    
    // 应用边界条件: 左固定，右拉伸
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 0.0});    // 左边界 u_x = 0
    bcs.push_back({"left", 1, 0.0});    // 左边界 u_y = 0
    bcs.push_back({"right", 0, 0.01});   // 右边界 u_x = 0.01 (拉伸!)
    
    FEM_INFO("\n=== Applying Boundary Conditions ===");
    FEM_INFO("BCs: left fixed, right pulled (u_x = 0.1)");
    
    // 检查受影响的DOF
    const auto& left_nodes = mesh.boundary("left");
    FEM_INFO("Left boundary nodes: " + std::to_string(left_nodes.size()));
    FEM_INFO("Constrained DOFs:");
    for (Index node_id : left_nodes) {
        Index dof_x = node_id * 2 + 0;
        Index dof_y = node_id * 2 + 1;
        FEM_INFO("  Node " + std::to_string(node_id) + " -> DOFs " + 
                 std::to_string(dof_x) + ", " + std::to_string(dof_y));
    }
    
    assembler.apply_dirichlet(bcs);
    
    // 获取边界条件后的矩阵
    SparseMatrixCSR K_after = assembler.matrix();
    const Vector& F_after = assembler.rhs();
    
    print_matrix_info(K_after, "Matrix After BC");
    
    FEM_INFO("\nRHS after BC (first 10):");
    for (std::size_t i = 0; i < std::min(F_after.size(), (std::size_t)10); ++i) {
        FEM_INFO("  F[" + std::to_string(i) + "] = " + fmt_sci(F_after[i]));
    }
    
    // 检查约束DOF的矩阵行
    FEM_INFO("\n=== Checking Constrained DOF Rows ===");
    Index check_dof = 0;  // DOF 0 是约束的
    FEM_INFO("Row " + std::to_string(check_dof) + " (constrained):");
    
    Index row_start = K_after.row_ptr()[check_dof];
    Index row_end = K_after.row_ptr()[check_dof + 1];
    
    for (Index k = row_start; k < row_end; ++k) {
        Index col = K_after.col_indices()[k];
        Real val = K_after.values()[k];
        
        if (std::abs(val) > 1e-12) {
            FEM_INFO("  K(" + std::to_string(check_dof) + "," + std::to_string(col) + ") = " + fmt_sci(val));
        }
    }
    
    // 检查一个未约束DOF的行
    FEM_INFO("\n=== Checking Free DOF Rows ===");
    Index free_dof = 2;  // DOF 2 应该是自由的
    FEM_INFO("Row " + std::to_string(free_dof) + " (free):");
    
    row_start = K_after.row_ptr()[free_dof];
    row_end = K_after.row_ptr()[free_dof + 1];
    
    int count = 0;
    for (Index k = row_start; k < row_end; ++k) {
        Index col = K_after.col_indices()[k];
        Real val = K_after.values()[k];
        
        if (std::abs(val) > 1e-12) {
            FEM_INFO("  K(" + std::to_string(free_dof) + "," + std::to_string(col) + ") = " + fmt_sci(val));
            count++;
            if (count > 10) {
                FEM_INFO("  ... (more entries)");
                break;
            }
        }
    }
    
    // 尝试求解
    FEM_INFO("\n=== Attempting to Solve ===");
    
    std::vector<Real> F_std = F_after.raw();
    std::vector<Real> u_std(F_after.size(), 0.0);
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(1000);
    
    auto result = solver.solve(K_after, F_std, u_std);
    
    FEM_INFO("Result: " + std::string(result.converged ? "CONVERGED" : "FAILED"));
    FEM_INFO("Residual: " + fmt_sci(result.residual));
    
    if (result.converged) {
        FEM_INFO("\n=== Solution ===");
        Real max_u = 0.0;
        for (std::size_t i = 0; i < u_std.size(); ++i) {
            max_u = std::max(max_u, std::abs(u_std[i]));
        }
        FEM_INFO("Max |u| = " + fmt_sci(max_u));
        FEM_INFO("Test PASSED!");
        return 0;
    } else {
        FEM_ERROR("Test FAILED!");
        return 1;
    }
}
