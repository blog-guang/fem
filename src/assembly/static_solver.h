/**
 * static_solver.h - Static Structural Solver
 * 
 * 高级静态结构求解器，集成：
 * - LoadStep 管理
 * - 增量求解
 * - 非线性求解（Newton-Raphson）
 * - 收敛历史记录
 */

#pragma once

#include "assembly/assembler.h"
#include "assembly/load_step.h"
#include "assembly/convergence_history.h"
#include "math/solver_factory.h"
#include "math/newton_raphson.h"
#include "physics/physics_base.h"
#include "mesh/model.h"
#include <memory>
#include <functional>

namespace fem {

/**
 * StaticSolver - 静态结构求解器
 * 
 * 用法示例：
 * ```cpp
 * StaticSolver solver(model, physics, dofs_per_node);
 * 
 * // 定义载荷步
 * LoadStep step1;
 * step1.set_displacement(node_id, DOF_Z, 0.0, 0.01);
 * step1.set_num_substeps(10);
 * solver.add_load_step(step1);
 * 
 * // 求解
 * solver.solve();
 * 
 * // 获取结果
 * Vector u = solver.solution();
 * solver.print_summary();
 * ```
 */
class StaticSolver {
public:
    // ═══ 构造 ═══
    
    /**
     * 构造静态求解器
     * 
     * @param model 有限元模型（网格+材料）
     * @param physics 物理模块
     * @param dofs_per_node 每节点自由度数
     */
    StaticSolver(const Model& model, 
                 physics::PhysicsBase* physics,
                 int dofs_per_node);
    
    ~StaticSolver() = default;
    
    // ═══ 载荷步管理 ═══
    
    /**
     * 添加载荷步
     */
    void add_load_step(const LoadStep& step);
    
    /**
     * 清空所有载荷步
     */
    void clear_load_steps();
    
    /**
     * 获取载荷步管理器
     */
    LoadStepManager& load_step_manager() { return load_steps_; }
    const LoadStepManager& load_step_manager() const { return load_steps_; }
    
    // ═══ 求解器设置 ═══
    
    /**
     * 设置线性求解器类型
     */
    void set_linear_solver(const std::string& solver_type,
                          const SolverFactory::Parameters& params = {},
                          const SolverFactory::StringParams& str_params = {});
    
    /**
     * 设置非线性求解器参数
     */
    void set_nonlinear_params(const NewtonRaphsonParams& params);
    
    /**
     * 启用/禁用非线性求解
     */
    void set_nonlinear(bool enable) { nonlinear_ = enable; }
    bool is_nonlinear() const { return nonlinear_; }
    
    // ═══ 求解 ═══
    
    /**
     * 执行求解（所有载荷步）
     * 
     * @return 是否全部收敛
     */
    bool solve();
    
    /**
     * 求解单个载荷步
     * 
     * @param load_step_id 载荷步 ID
     * @return 是否收敛
     */
    bool solve_load_step(int load_step_id);
    
    /**
     * 求解单个子步（低级接口）
     * 
     * @param step 载荷步
     * @param substep_id 子步 ID
     * @return 是否收敛
     */
    bool solve_substep(const LoadStep& step, int substep_id);
    
    // ═══ 结果访问 ═══
    
    /**
     * 获取当前解向量
     */
    const Vector& solution() const { return u_; }
    Vector& solution() { return u_; }
    
    /**
     * 获取收敛历史
     */
    const ConvergenceHistory& convergence_history() const { return history_; }
    ConvergenceHistory& convergence_history() { return history_; }
    
    /**
     * 获取装配器（用于后处理）
     */
    const Assembler& assembler() const { return assembler_; }
    
    // ═══ 工具函数 ═══
    
    /**
     * 打印求解摘要
     */
    void print_summary() const;
    
    /**
     * 导出收敛历史到 CSV
     */
    void export_history(const std::string& filename) const;
    
private:
    // ═══ 数据成员 ═══
    const Model& model_;
    physics::PhysicsBase* physics_;
    int dofs_per_node_;
    
    Assembler assembler_;
    LoadStepManager load_steps_;
    ConvergenceHistory history_;
    
    Vector u_;                      // 当前解向量
    bool nonlinear_;                // 是否非线性求解
    
    // 求解器
    std::unique_ptr<LinearSolver> linear_solver_;
    NewtonRaphsonParams nr_params_;
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 线性求解（单个子步）
     */
    bool solve_linear(const std::map<std::pair<Index, int>, Real>& bcs,
                     const std::map<std::pair<Index, int>, Real>& forces);
    
    /**
     * 非线性求解（单个子步）
     */
    bool solve_nonlinear(const std::map<std::pair<Index, int>, Real>& bcs,
                        const std::map<std::pair<Index, int>, Real>& forces);
};

}  // namespace fem
