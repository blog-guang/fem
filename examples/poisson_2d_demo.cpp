/**
 * poisson_2d_demo.cpp
 * 
 * 2D Poisson 方程 -∇·(k∇u) = f 的有限元求解演示
 * 
 * 问题设定:
 * - 区域: [0,1]x[0,1]
 * - u = 0 (Dirichlet 边界条件)
 * - f = 1 (源项)
 * - k = 1 (导热系数)
 * 
 * 解析解: u(x,y) = (x(1-x)y(1-y))/8 (近似)
 * 
 * 流程: 网格生成 → 单元装配 → 边界条件 → 求解
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "math/sparse_matrix.h"
#include "assembly/sparse_matrix.h"
#include "solver/cg.h"
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
    grad[0][0] = (p1[1] - p2[1]) / (2.0 * area);  // ∂N0/∂x
    grad[0][1] = (p2[0] - p1[0]) / (2.0 * area);  // ∂N0/∂y
    
    // ∇N1
    grad[1][0] = (p2[1] - p0[1]) / (2.0 * area);
    grad[1][1] = (p0[0] - p2[0]) / (2.0 * area);
    
    // ∇N2
    grad[2][0] = (p0[1] - p1[1]) / (2.0 * area);
    grad[2][1] = (p1[0] - p0[0]) / (2.0 * area);
}

int main() {
    FEM_INFO("=== 2D Poisson Equation Demo ===");
    
    // ── 1. 创建模型 ──
    Model model("Poisson 2D");
    
    // 材料 (导热系数)
    int mat_id = model.add_material("Conductivity");
    model.material(mat_id).set_property("k", 1.0);  // 导热系数
    
    // 网格 (2D 三角形网格)
    int mesh_id = model.add_mesh("domain", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 20x20 网格
    MeshGenerator::generate_unit_square_tri(20, 20, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh generated: " + std::to_string(mesh.num_nodes()) + 
             " nodes, " + std::to_string(mesh.num_elements()) + " elements");
    
    // ── 2. 装配刚度矩阵和载荷向量 ──
    Timer timer;
    timer.start();
    
    std::size_t n_nodes = mesh.num_nodes();
    
    // 使用 COO 格式装配 (适合装配阶段)
    SparseMatrixCOO K_coo(n_nodes, n_nodes);
    Vector F(n_nodes, 0.0);
    
    // 遌溉刚度矩阵和载荷向量
    for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        const Element& elem = mesh.element(elem_id);
        
        if (elem.type() != ElementType::Tri3) continue;  // 只处理 Tri3
        
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
        
        DenseMatrix Ke(3, 3);
        Real k = mesh.material()->property("k", 1.0);  // 导热系数
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // Ke[i][j] = k * ∫_Ω ∇Ni · ∇Nj dΩ = k * (∇Ni · ∇Nj) * area
                Ke(i, j) = k * (grad[i][0]*grad[j][0] + grad[i][1]*grad[j][1]) * area;
            }
        }
        
        // 计算单元载荷向量 Fe[i] = ∫_Ω f * Ni dΩ (f=1)
        Real source = 1.0;  // 源项
        Vector Fe(3);  // 单元载荷向量
        for (int i = 0; i < 3; ++i) {
            Fe[i] = source * area / 3.0;  // 每个节点承担 1/3 的面积载荷
        }
        
        // 添加到全局矩阵和向量
        for (int i = 0; i < 3; ++i) {
            Index gi = nodes[i];  // 全局节点编号
            
            F[gi] += Fe[i];
            
            for (int j = 0; j < 3; ++j) {
                Index gj = nodes[j];  // 全局节点编号
                K_coo.add(gi, gj, Ke(i, j));
            }
        }
    }
    
    FEM_INFO("Assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 3. 转换为 CSR 格式 (适合求解) ──
    timer.start();
    SparseMatrixCSR K_new = coo_to_csr(K_coo);
    
    // 转换到旧格式供求解器使用
    CSRMatrix K;
    K.rows = K_new.rows();
    K.row_ptr = K_new.row_ptr();
    K.col_idx = K_new.col_indices();
    K.values = K_new.values();
    
    FEM_INFO("Matrix format conversion in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 4. 施加边界条件 (u = 0 on all boundaries) ──
    timer.start();
    
    // 获取边界节点
    std::vector<bool> is_boundary(n_nodes, false);
    for (const std::string& name : {"left", "right", "bottom", "top"}) {
        const auto& boundary_nodes = mesh.boundary(name);
        for (Index node_id : boundary_nodes) {
            if (node_id < n_nodes) {
                is_boundary[node_id] = true;
            }
        }
    }
    
    // 对边界节点施加 Dirichlet 条件 u = 0
    for (std::size_t i = 0; i < n_nodes; ++i) {
        if (is_boundary[i]) {
            // 设置 K(ii) = 1, K(i,j) = 0 (j≠i), F(i) = 0
            // 找到第i行的起始位置
            Index row_start = K.row_ptr[i];
            Index row_end = K.row_ptr[i + 1];
            
            for (Index k = row_start; k < row_end; ++k) {
                Index j = K.col_idx[k];
                if (j == i) {
                    K.values[k] = 1.0;
                } else {
                    K.values[k] = 0.0;
                }
            }
            F[i] = 0.0;
        }
    }
    
    FEM_INFO("Boundary conditions applied in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 5. 求解线性系统 Ku = F ──
    timer.start();
    
    Vector u(n_nodes, 0.0);  // 解向量
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(1000);
    
    // 转换 Vector 到 std::vector<Real>
    std::vector<Real> F_std = F.raw();
    std::vector<Real> u_std = u.raw();
    
    auto result = solver.solve(K, F_std, u_std);
    
    // 将结果转回 fem::Vector
    for (std::size_t i = 0; i < u.size(); ++i) {
        u[i] = u_std[i];
    }
    
    FEM_INFO("Solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("CG convergence: " + std::string(result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(result.residual));
    
    if (result.converged) {
        // ── 6. 输出结果 ──
        Real max_u = 0.0;
        std::size_t max_idx = 0;
        for (std::size_t i = 0; i < n_nodes; ++i) {
            if (std::abs(u[i]) > max_u) {
                max_u = std::abs(u[i]);
                max_idx = i;
            }
        }
        
        const auto& max_coord = mesh.node(max_idx).coords();
        FEM_INFO("Max solution: u = " + fmt_sci(max_u) + 
                 " at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
        
        // 与解析解比较 (中心点附近)
        Real center_exact = 1.0 / 64.0;  // 近似值
        FEM_INFO("Center approx: " + fmt_sci(center_exact));
        
        FEM_INFO("Demo completed successfully!");
    } else {
        FEM_ERROR("Solver failed to converge");
        return 1;
    }
    
    return 0;
}