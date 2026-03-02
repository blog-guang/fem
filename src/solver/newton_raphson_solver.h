/**
 * newton_raphson_solver.h - Newton-Raphson 非线性求解器
 * 
 * 用于几何非线性、材料非线性问题
 * 
 * 算法：
 * 1. 初始化：u⁰ = 0
 * 2. 迭代：
 *    - 计算残差：R = F_int(u^k) - F_ext
 *    - 计算切线刚度：K_t = ∂R/∂u
 *    - 求解：K_t Δu = -R
 *    - 更新：u^(k+1) = u^k + Δu
 *    - 检查收敛：||R|| < tol 或 ||Δu|| < tol
 */

#pragma once

#include "core/types.h"
#include "math/solver.h"
#include "math/sparse_matrix.h"
#include "assembly/assembler.h"
#include "physics/physics_base.h"
#include "mesh/model.h"
#include <functional>
#include <memory>

namespace fem {
namespace solver {

/**
 * Newton-Raphson 收敛准则
 */
struct NewtonRaphsonCriteria {
    Real residual_tol = 1e-6;       ///< 残差范数容差
    Real displacement_tol = 1e-8;   ///< 位移增量容差
    Real energy_tol = 1e-10;        ///< 能量容差
    int max_iterations = 50;        ///< 最大迭代次数
    
    bool use_residual = true;       ///< 启用残差判据
    bool use_displacement = true;   ///< 启用位移判据
    bool use_energy = false;        ///< 启用能量判据
};

/**
 * Newton-Raphson 迭代信息
 */
struct NewtonRaphsonIteration {
    int iteration = 0;              ///< 当前迭代次数
    Real residual_norm = 0.0;       ///< 残差范数
    Real displacement_norm = 0.0;   ///< 位移增量范数
    Real energy = 0.0;              ///< 能量 ΔuᵀR
    bool converged = false;         ///< 是否收敛
};

/**
 * Newton-Raphson 非线性求解器
 * 
 * 求解非线性方程：R(u) = F_int(u) - F_ext = 0
 */
class NewtonRaphsonSolver {
public:
    /**
     * 构造函数
     * 
     * @param model 有限元模型
     * @param physics 物理模块（必须支持几何非线性）
     * @param dofs_per_node 每节点自由度数
     */
    NewtonRaphsonSolver(
        Model& model,
        physics::PhysicsBase* physics,
        int dofs_per_node);
    
    // ═══ 求解控制 ═══
    
    /**
     * 求解非线性问题
     * 
     * @param F_ext 外力向量
     * @param u_initial 初始位移（默认零向量）
     * @return 是否收敛
     */
    bool solve(const Vector& F_ext, const Vector& u_initial = Vector());
    
    /**
     * 单步 Newton-Raphson 迭代
     * 
     * @return 是否收敛
     */
    bool iterate();
    
    // ═══ 结果访问 ═══
    
    /**
     * 获取当前解
     */
    const Vector& solution() const { return u_current_; }
    
    /**
     * 获取迭代历史
     */
    const std::vector<NewtonRaphsonIteration>& iteration_history() const {
        return history_;
    }
    
    /**
     * 获取最后一次迭代信息
     */
    const NewtonRaphsonIteration& last_iteration() const {
        return history_.back();
    }
    
    // ═══ 配置 ═══
    
    /**
     * 设置收敛准则
     */
    void set_criteria(const NewtonRaphsonCriteria& criteria) {
        criteria_ = criteria;
    }
    
    /**
     * 获取收敛准则
     */
    const NewtonRaphsonCriteria& criteria() const { return criteria_; }
    
    /**
     * 设置线性求解器
     */
    void set_linear_solver(std::unique_ptr<LinearSolver> solver) {
        linear_solver_ = std::move(solver);
    }
    
    /**
     * 启用/禁用详细输出
     */
    void set_verbose(bool verbose) { verbose_ = verbose; }
    
    // ═══ 高级功能 ═══
    
    /**
     * 设置载荷缩放因子（增量加载）
     * 
     * @param scale 载荷缩放因子 [0, 1]
     */
    void set_load_scale(Real scale) { load_scale_ = scale; }
    
    /**
     * 获取当前载荷缩放因子
     */
    Real load_scale() const { return load_scale_; }
    
    /**
     * 重置求解器（清除历史）
     */
    void reset();
    
private:
    // ═══ 内部计算 ═══
    
    /**
     * 计算残差：R = F_int - F_ext
     */
    void compute_residual(const Vector& u, Vector& R);
    
    /**
     * 计算内力向量
     */
    void compute_internal_force(const Vector& u, Vector& F_int);
    
    /**
     * 计算切线刚度矩阵
     */
    void compute_tangent_stiffness(const Vector& u, SparseMatrixCSR& K_t);
    
    /**
     * 应用边界条件到残差和切线刚度
     */
    void apply_boundary_conditions(Vector& R, SparseMatrixCSR& K_t);
    
    /**
     * 检查收敛
     */
    bool check_convergence(
        const Vector& R,
        const Vector& du,
        NewtonRaphsonIteration& iter);
    
    /**
     * 打印迭代信息
     */
    void print_iteration(const NewtonRaphsonIteration& iter);
    
    // ═══ 成员变量 ═══
    
    Model& model_;                              ///< 有限元模型
    physics::PhysicsBase* physics_;             ///< 物理模块
    int dofs_per_node_;                         ///< 每节点自由度数
    
    NewtonRaphsonCriteria criteria_;            ///< 收敛准则
    std::unique_ptr<LinearSolver> linear_solver_; ///< 线性求解器
    bool verbose_;                              ///< 详细输出
    
    Vector u_current_;                          ///< 当前位移
    Vector F_ext_;                              ///< 外力向量
    Real load_scale_;                           ///< 载荷缩放因子
    
    std::vector<NewtonRaphsonIteration> history_; ///< 迭代历史
    
    std::unique_ptr<Assembler> assembler_;      ///< 装配器
};

}  // namespace solver
}  // namespace fem
