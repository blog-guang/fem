/**
 * nonlinear_truss.cpp - 几何非线性桁架示例
 * 
 * 问题：两杆桁架在大变形下的非线性响应
 * 
 * 几何：
 *     节点 2 (自由)
 *       /  \
 *      /    \
 *   节点0   节点1
 *  (固定)  (固定)
 * 
 * 载荷：节点2受垂直向下力 P
 * 
 * 非线性来源：
 * - 几何非线性：大位移导致杆件方向改变
 * - 应变-位移关系非线性：ε = (L' - L)/L
 * 
 * 材料：线弹性 (σ = E*ε)
 */

#include "solver/newton_raphson.h"
#include "solver/cg.h"
#include "math/sparse_matrix.h"
#include "core/logger.h"
#include "core/timer.h"

#include <cmath>
#include <vector>

using namespace fem;

/**
 * 两杆桁架几何非线性问题
 * 
 * 自由度：节点2的位移 (u_x, u_y)
 * 残差：R = F_int - F_ext
 */
class NonlinearTruss : public NonlinearProblem {
public:
    NonlinearTruss(Real E, Real A, Real L, Real angle, Real P)
        : E_(E), A_(A), L_(L), angle_(angle), P_(P) {
        
        // 初始节点坐标
        node0_x_ = 0.0;
        node0_y_ = 0.0;
        node1_x_ = 2.0 * L * std::cos(angle);
        node1_y_ = 0.0;
        node2_x0_ = L * std::cos(angle);
        node2_y0_ = L * std::sin(angle);
    }
    
    void compute_residual(const std::vector<Real>& u, 
                         std::vector<Real>& R) override {
        R.resize(2);
        
        // 当前节点2位置
        Real x2 = node2_x0_ + u[0];
        Real y2 = node2_y0_ + u[1];
        
        // 计算两根杆的内力
        Real F0x, F0y, F1x, F1y;
        compute_bar_force(node0_x_, node0_y_, x2, y2, F0x, F0y);
        compute_bar_force(node1_x_, node1_y_, x2, y2, F1x, F1y);
        
        // 残差 = F_int - F_ext
        R[0] = F0x + F1x;          // x方向内力平衡
        R[1] = F0y + F1y - P_;     // y方向内力平衡 (减去外力)
    }
    
    void compute_tangent(const std::vector<Real>& u,
                        SparseMatrixCSR& K_t) override {
        // 当前节点2位置
        Real x2 = node2_x0_ + u[0];
        Real y2 = node2_y0_ + u[1];
        
        // 计算两根杆的切线刚度贡献
        Real K0[2][2], K1[2][2];
        compute_bar_tangent(node0_x_, node0_y_, x2, y2, K0);
        compute_bar_tangent(node1_x_, node1_y_, x2, y2, K1);
        
        // 组装总切线刚度
        SparseMatrixCOO K_coo(2, 2);
        K_coo.add(0, 0, K0[0][0] + K1[0][0]);
        K_coo.add(0, 1, K0[0][1] + K1[0][1]);
        K_coo.add(1, 0, K0[1][0] + K1[1][0]);
        K_coo.add(1, 1, K0[1][1] + K1[1][1]);
        
        K_t = coo_to_csr(K_coo);
    }
    
private:
    Real E_;     // 杨氏模量
    Real A_;     // 截面积
    Real L_;     // 初始长度
    Real angle_; // 初始角度
    Real P_;     // 外力
    
    // 节点初始坐标
    Real node0_x_, node0_y_;
    Real node1_x_, node1_y_;
    Real node2_x0_, node2_y0_;
    
    /**
     * 计算杆件内力
     * 
     * @param x1, y1 杆件起点
     * @param x2, y2 杆件终点
     * @param Fx, Fy 作用在终点的内力
     */
    void compute_bar_force(Real x1, Real y1, Real x2, Real y2,
                          Real& Fx, Real& Fy) {
        // 当前杆长和方向
        Real dx = x2 - x1;
        Real dy = y2 - y1;
        Real L_current = std::sqrt(dx*dx + dy*dy);
        
        // 应变（工程应变）
        Real strain = (L_current - L_) / L_;
        
        // 应力
        Real stress = E_ * strain;
        
        // 轴力
        Real force = stress * A_;
        
        // 内力分量（指向杆内）
        Fx = force * dx / L_current;
        Fy = force * dy / L_current;
    }
    
    /**
     * 计算杆件切线刚度矩阵 (2x2)
     * 
     * K = (EA/L) * (1 + ε) * (n⊗n) + (EA/L²) * ε * (I - n⊗n)
     * 
     * 其中：
     * - n: 杆件方向向量
     * - ε: 应变
     * - ⊗: 外积
     */
    void compute_bar_tangent(Real x1, Real y1, Real x2, Real y2,
                            Real K[2][2]) {
        // 当前杆长和方向
        Real dx = x2 - x1;
        Real dy = y2 - y1;
        Real L_current = std::sqrt(dx*dx + dy*dy);
        
        // 方向余弦
        Real nx = dx / L_current;
        Real ny = dy / L_current;
        
        // 应变
        Real strain = (L_current - L_) / L_;
        
        // 几何刚度和材料刚度
        Real k_material = E_ * A_ * (1.0 + strain) / L_current;
        Real k_geometric = E_ * A_ * strain / (L_current * L_current);
        
        // 切线刚度矩阵
        K[0][0] = k_material * nx * nx + k_geometric * (1.0 - nx * nx);
        K[0][1] = k_material * nx * ny - k_geometric * nx * ny;
        K[1][0] = k_material * ny * nx - k_geometric * ny * nx;
        K[1][1] = k_material * ny * ny + k_geometric * (1.0 - ny * ny);
    }
};

int main() {
    FEM_INFO("=== Geometrically Nonlinear Truss Demo ===");
    
    // ── 问题参数 ──
    Real E = 2e5;           // 杨氏模量 (MPa)
    Real A = 100.0;         // 截面积 (mm²)
    Real L = 1000.0;        // 初始杆长 (mm)
    Real angle = M_PI / 6;  // 初始角度 (30°)
    
    // 载荷步
    std::vector<Real> load_steps = {10.0, 50.0, 100.0, 200.0, 500.0};
    
    FEM_INFO("Truss parameters:");
    FEM_INFO("  E = " + fmt_sci(E) + " MPa");
    FEM_INFO("  A = " + std::to_string(A) + " mm²");
    FEM_INFO("  L = " + std::to_string(L) + " mm");
    FEM_INFO("  angle = " + std::to_string(angle * 180.0 / M_PI) + "°");
    
    // ── 载荷递增求解 ──
    std::vector<Real> u = {0.0, 0.0};  // 累积位移
    
    Timer timer;
    
    for (Real P : load_steps) {
        FEM_INFO("\n=== Load Step: P = " + std::to_string(P) + " N ===");
        timer.start();
        
        // 创建非线性问题
        NonlinearTruss problem(E, A, L, angle, P);
        
        // Newton-Raphson 求解
        NewtonRaphsonSolver solver;
        NewtonRaphsonParams params;
        params.tol = 1e-6;
        params.tol_relative = 1e-8;
        params.max_iter = 20;
        params.verbose = true;
        solver.set_params(params);
        
        auto result = solver.solve(problem, u);
        
        Real solve_time = timer.elapsed_s();
        
        if (result.converged) {
            FEM_INFO("Converged in " + std::to_string(result.iterations) + 
                    " iterations, " + std::to_string(solve_time * 1000.0) + " ms");
            FEM_INFO("Solution:");
            FEM_INFO("  u_x = " + fmt_sci(u[0]) + " mm");
            FEM_INFO("  u_y = " + fmt_sci(u[1]) + " mm");
            FEM_INFO("  |u| = " + fmt_sci(std::sqrt(u[0]*u[0] + u[1]*u[1])) + " mm");
        } else {
            FEM_ERROR("Failed to converge at load P = " + std::to_string(P));
            break;
        }
    }
    
    FEM_INFO("\n=== Nonlinear Analysis Complete ===");
    FEM_INFO("Geometric nonlinearity successfully captured!");
    FEM_INFO("Large displacement: " + fmt_sci(std::sqrt(u[0]*u[0] + u[1]*u[1])) + " mm");
    
    return 0;
}
