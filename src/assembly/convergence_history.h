/**
 * convergence_history.h - Convergence History Tracking
 * 
 * 记录非线性求解的收敛历史
 */

#pragma once

#include "core/types.h"
#include <vector>
#include <string>

namespace fem {

/**
 * SubstepResult - 单个子步的求解结果
 */
struct SubstepResult {
    int load_step_id;           // 载荷步 ID
    int substep_id;             // 子步 ID
    Real time;                  // 伪时间
    bool converged;             // 是否收敛
    int iterations;             // 迭代次数
    Real residual_norm;         // 残差范数
    Real displacement_norm;     // 位移范数
    Real solve_time;            // 求解时间（秒）
};

/**
 * ConvergenceHistory - 收敛历史记录
 * 
 * 记录所有载荷步和子步的收敛信息
 */
class ConvergenceHistory {
public:
    /**
     * 添加子步结果
     */
    void add_substep(const SubstepResult& result);
    
    /**
     * 获取所有子步结果
     */
    const std::vector<SubstepResult>& substeps() const { return substeps_; }
    
    /**
     * 获取总迭代次数
     */
    int total_iterations() const;
    
    /**
     * 获取总求解时间
     */
    Real total_time() const;
    
    /**
     * 是否全部收敛
     */
    bool all_converged() const;
    
    /**
     * 打印摘要
     */
    void print_summary() const;
    
    /**
     * 导出到文件（CSV 格式）
     */
    void export_csv(const std::string& filename) const;
    
    /**
     * 清空历史
     */
    void clear() { substeps_.clear(); }
    
private:
    std::vector<SubstepResult> substeps_;
};

}  // namespace fem
