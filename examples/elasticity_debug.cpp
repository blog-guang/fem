/**
 * elasticity_debug.cpp - Elasticity2D 详细调试
 * 
 * 单个三角形单元测试
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_v2.h"
#include "math/dense_matrix.h"
#include "mesh/element.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== Elasticity2D Debug ===");
    
    // 创建最简单的模型：1个三角形
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 手动添加3个节点和1个单元
    mesh.add_node({0.0, 0.0, 0.0});
    mesh.add_node({1.0, 0.0, 0.0});
    mesh.add_node({0.0, 1.0, 0.0});
    mesh.add_element<Tri3>(0, 1, 2);
    
    FEM_INFO("Single element mesh: 3 nodes, 1 element");
    
    // 创建物理模块
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    // 手动计算单元矩阵
    DenseMatrix Ke(6, 6);
    Vector Fe(6);
    
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 手动计算 B 矩阵和 D 矩阵来验证
    FEM_INFO("\n=== Manual B Matrix Calculation ===");
    Vec3 coords[3];
    coords[0] = mesh.node(0).coords();
    coords[1] = mesh.node(1).coords();
    coords[2] = mesh.node(2).coords();
    
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];
    
    Real detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    Real area = 0.5 * std::abs(detJ);
    
    FEM_INFO("Nodes: (" + fmt_sci(x0) + "," + fmt_sci(y0) + "), (" + 
             fmt_sci(x1) + "," + fmt_sci(y1) + "), (" + 
             fmt_sci(x2) + "," + fmt_sci(y2) + ")");
    FEM_INFO("Area: " + fmt_sci(area));
    
    Real dN_dx[3], dN_dy[3];
    dN_dx[0] = (y1 - y2) / (2.0 * area);
    dN_dy[0] = (x2 - x1) / (2.0 * area);
    dN_dx[1] = (y2 - y0) / (2.0 * area);
    dN_dy[1] = (x0 - x2) / (2.0 * area);
    dN_dx[2] = (y0 - y1) / (2.0 * area);
    dN_dy[2] = (x1 - x0) / (2.0 * area);
    
    FEM_INFO("Gradients:");
    for (int i = 0; i < 3; ++i) {
        FEM_INFO("  dN" + std::to_string(i) + "/dx = " + fmt_sci(dN_dx[i]) + 
                 ", dN" + std::to_string(i) + "/dy = " + fmt_sci(dN_dy[i]));
    }
    
    DenseMatrix B(3, 6, 0.0);
    for (int i = 0; i < 3; ++i) {
        int col_u = i * 2;
        int col_v = i * 2 + 1;
        B(0, col_u) = dN_dx[i];
        B(1, col_v) = dN_dy[i];
        B(2, col_u) = dN_dy[i];
        B(2, col_v) = dN_dx[i];
    }
    
    FEM_INFO("\nB Matrix (3x6):");
    for (int i = 0; i < 3; ++i) {
        std::string row = "Row " + std::to_string(i) + ": ";
        for (int j = 0; j < 6; ++j) {
            char buf[16];
            std::snprintf(buf, sizeof(buf), "%8.3f ", B(i, j));
            row += buf;
        }
        FEM_INFO(row);
    }
    
    // D 矩阵 (平面应力)
    Real E = 1000.0;
    Real nu = 0.3;
    Real factor = E / (1.0 - nu * nu);
    DenseMatrix D(3, 3, 0.0);
    D(0, 0) = factor;
    D(0, 1) = factor * nu;
    D(1, 0) = factor * nu;
    D(1, 1) = factor;
    D(2, 2) = factor * (1.0 - nu) / 2.0;
    
    FEM_INFO("\nD Matrix (3x3) - Plane Stress:");
    for (int i = 0; i < 3; ++i) {
        std::string row = "Row " + std::to_string(i) + ": ";
        for (int j = 0; j < 3; ++j) {
            char buf[16];
            std::snprintf(buf, sizeof(buf), "%10.2f ", D(i, j));
            row += buf;
        }
        FEM_INFO(row);
    }
    
    // 手动计算 Ke = B^T D B * area
    DenseMatrix DB = D * B;  // 3x6
    DenseMatrix Bt = B.transpose();  // 6x3
    DenseMatrix Ke_manual = Bt * DB * area;  // 6x6
    
    FEM_INFO("\nManual Ke (6x6):");
    for (int i = 0; i < 6; ++i) {
        std::string row = "Row " + std::to_string(i) + ": ";
        for (int j = 0; j < 6; ++j) {
            char buf[16];
            std::snprintf(buf, sizeof(buf), "%10.2f ", Ke_manual(i, j));
            row += buf;
        }
        FEM_INFO(row);
    }
    
    FEM_INFO("=== Element Matrix Ke (6x6) ===");
    for (int i = 0; i < 6; ++i) {
        std::string row = "Row " + std::to_string(i) + ": ";
        for (int j = 0; j < 6; ++j) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%10.3e ", Ke(i, j));
            row += buf;
        }
        FEM_INFO(row);
    }
    
    // 检查对称性
    FEM_INFO("\n=== Symmetry Check ===");
    Real max_asym = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = i+1; j < 6; ++j) {
            Real diff = std::abs(Ke(i, j) - Ke(j, i));
            max_asym = std::max(max_asym, diff);
            if (diff > 1e-10) {
                FEM_WARN("Asymmetry at (" + std::to_string(i) + "," + std::to_string(j) + "): " + 
                         std::to_string(diff));
            }
        }
    }
    FEM_INFO("Max asymmetry: " + fmt_sci(max_asym));
    
    // 检查正定性
    FEM_INFO("\n=== Positive Definiteness Check ===");
    bool all_positive = true;
    for (int i = 0; i < 6; ++i) {
        if (Ke(i, i) <= 0.0) {
            FEM_ERROR("Diagonal element K(" + std::to_string(i) + "," + std::to_string(i) + 
                     ") = " + std::to_string(Ke(i, i)) + " <= 0");
            all_positive = false;
        }
    }
    
    if (all_positive) {
        FEM_INFO("All diagonal elements are positive ✓");
    }
    
    // 检查是否有 NaN/Inf
    FEM_INFO("\n=== NaN/Inf Check ===");
    bool has_nan = false;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (std::isnan(Ke(i, j)) || std::isinf(Ke(i, j))) {
                FEM_ERROR("K(" + std::to_string(i) + "," + std::to_string(j) + ") = " + 
                         std::to_string(Ke(i, j)) + " is NaN/Inf");
                has_nan = true;
            }
        }
    }
    
    if (!has_nan) {
        FEM_INFO("No NaN/Inf in element matrix ✓");
    } else {
        FEM_ERROR("Element matrix contains NaN/Inf!");
        return 1;
    }
    
    // 测试 v^T K v > 0 (正定性)
    FEM_INFO("\n=== Positive Definite Test (v^T K v) ===");
    Vector v(6, 1.0);  // 全1向量
    
    Real vKv = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            vKv += v[i] * Ke(i, j) * v[j];
        }
    }
    
    FEM_INFO("v^T K v = " + fmt_sci(vKv));
    
    if (vKv > 0.0) {
        FEM_INFO("Positive definite ✓");
    } else {
        FEM_ERROR("NOT positive definite!");
        return 1;
    }
    
    FEM_INFO("\n=== Single Element Test PASSED ===");
    return 0;
}
