/**
 * nonlinear_solver.h - 非线性求解器（使用回调函数）
 * 
 * 集成 Newton-Raphson 方法与现有 Assembler API
 */

#pragma once

#include "core/types.h"
#include "assembly/assembler.h"
#include "physics/physics_base.h"
#include "physics/elasticity_nonlinear.h"
#include "mesh/model.h"
#include "math/solver.h"
#include <functional>
#include <memory>
#include <vector>

namespace fem {
namespace solver {

/**
 * 非线性求解器收敛准则
 */
struct NonlinearCriteria {
    Real residual_tol = 1e-6;       ///< 残差范数容差
    Real displacement_tol = 1e-8;   ///< 位移增量容差
    Real energy_tol = 1e-10;        ///< 能量容差（ΔuᵀR）
    int max_iterations = 50;        ///< 最大迭代次数
    
    bool use_residual = true;       ///< 启用残差判据
    bool use_displacement = true;   ///< 启用位移判据
    bool use_energy = false;        ///< 启用能量判据
};

/**
 * 非线性迭代信息
 */
struct NonlinearIteration {
    int iteration = 0;              ///< 当前迭代次数
    Real residual_norm = 0.0;       ///< 残差范数 ||R||
    Real displacement_norm = 0.0;   ///< 位移增量范数 ||Δu||
    Real energy = 0.0;              ///< 能量 |ΔuᵀR|
    bool converged = false;         ///< 是否收敛
    Real solve_time = 0.0;          ///< 线性求解时间（秒）
};

/**
 * 非线性求解器（Newton-Raphson 方法）
 * 
 * 使用回调函数与 Assembler 集成
 * 
 * 用法：
 * ```cpp
 * NonlinearSolver solver(model);
 * solver.set_physics(&physics);
 * solver.set_criteria(criteria);
 * 
 * bool converged = solver.solve(F_ext);
 * Vector u = solver.solution();
 * ```
 */
class NonlinearSolver {
public:
    /**
     * 构造函数
     * 
     * @param model 有限元模型
     * @param dofs_per_node 每节点自由度数
     */
    NonlinearSolver(const Model& model, Index dofs_per_node = 3);
    
    // ═══ 配置 ═══
    
    /**
     * 设置物理模块（必须支持几何非线性）
     */
    void set_physics(physics::ElasticityNonlinear* physics) {
        physics_ = physics;
    }
    
    /**
     * 设置收敛准则
     */
    void set_criteria(const NonlinearCriteria& criteria) {
        criteria_ = criteria;
    }
    
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
    
    /**
     * 设置载荷缩放因子（用于增量加载）
     */
    void set_load_scale(Real scale) { load_scale_ = scale; }
    
    // ═══ 边界条件 ═══
    
    /**
     * 应用 Dirichlet 边界条件
     * 
     * @param node_id 节点 ID
     * @param dof 自由度索引（0=x, 1=y, 2=z）
     * @param value 约束值
     */
    void apply_dirichlet(Index node_id, Index dof, Real value);
    
    /**
     * 批量应用 Dirichlet 边界条件
     */
    void apply_dirichlet(const std::vector<DirichletBC>& bcs);
    
    // ═══ 求解 ═══
    
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
    
    /**
     * 重置求解器（清除历史和边界条件）
     */
    void reset();
    
    // ═══ 结果访问 ═══
    
    /**
     * 获取当前解
     */
    const Vector& solution() const { return u_current_; }
    
    /**
     * 获取迭代历史
     */
    const std::vector<NonlinearIteration>& history() const {
        return history_;
    }
    
    /**
     * 获取最后一次迭代信息
     */
    const NonlinearIteration& last_iteration() const {
        return history_.back();
    }
    
    /**
     * 获取迭代次数
     */
    int num_iterations() const {
        return static_cast<int>(history_.size());
    }
    
private:
    // ═══ 内部方法 ═══
    
    /**
     * 计算残差：R = F_int - F_ext
     */
    void compute_residual(const Vector& u, Vector& R);
    
    /**
     * 计算内力向量（使用物理模块）
     */
    void compute_internal_force(const Vector& u, Vector& F_int);
    
    /**
     * 计算切线刚度矩阵（使用 Assembler + 回调）
     */
    void compute_tangent_stiffness(const Vector& u, SparseMatrixCSR& K_t);
    
    /**
     * 应用边界条件到刚度矩阵和残差向量
     */
    void apply_bc_to_system(SparseMatrixCSR& K, Vector& R);
    
    /**
     * 检查收敛
     */
    bool check_convergence(
        const Vector& R,
        const Vector& du,
        NonlinearIteration& iter);
    
    /**
     * 打印迭代信息
     */
    void print_iteration(const NonlinearIteration& iter);
    
    // ═══ 成员变量 ═══
    
    const Model& model_;                        ///< 有限元模型
    Index dofs_per_node_;                       ///< 每节点自由度数
    Index total_dofs_;                          ///< 总自由度数
    
    physics::ElasticityNonlinear* physics_;     ///< 物理模块
    NonlinearCriteria criteria_;                ///< 收敛准则
    std::unique_ptr<LinearSolver> linear_solver_; ///< 线性求解器
    bool verbose_;                              ///< 详细输出
    Real load_scale_;                           ///< 载荷缩放因子
    
    Vector u_current_;                          ///< 当前位移
    Vector F_ext_;                              ///< 外力向量
    std::vector<NonlinearIteration> history_;   ///< 迭代历史
    
    // Dirichlet 边界条件存储
    std::map<Index, Real> dirichlet_bc_;        ///< (global_dof, value)
};

}  // namespace solver
}  // namespace fem
