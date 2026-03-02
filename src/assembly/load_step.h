/**
 * load_step.h - Load Step and Substep Management
 * 
 * 载荷步管理：支持增量求解、载荷递增/阶跃
 */

#pragma once

#include "core/types.h"
#include <map>
#include <vector>
#include <string>

namespace fem {

/**
 * 载荷施加方式
 */
enum class LoadRampType {
    RAMPED,   // 线性递增（默认）
    STEPPED   // 阶跃（瞬间施加）
};

/**
 * LoadStep - 载荷步
 * 
 * 一个完整的载荷历史，可以细分为多个子步（substeps）
 * 
 * 用法示例：
 * ```cpp
 * LoadStep step1;
 * step1.set_displacement(node_id, DOF_X, 0.0, 0.01);  // 0 → 0.01 m
 * step1.set_num_substeps(10);
 * step1.set_ramp_type(LoadRampType::RAMPED);
 * ```
 */
class LoadStep {
public:
    // ═══ 构造 ═══
    LoadStep(int id = 1) 
        : id_(id), num_substeps_(1), ramp_type_(LoadRampType::RAMPED),
          time_start_(0.0), time_end_(1.0), load_scale_(1.0) {}
    
    // ═══ 基本设置 ═══
    
    /**
     * 设置载荷步 ID
     */
    void set_id(int id) { id_ = id; }
    int id() const { return id_; }
    
    /**
     * 设置子步数量
     */
    void set_num_substeps(int num) { 
        if (num < 1) num = 1;
        num_substeps_ = num; 
    }
    int num_substeps() const { return num_substeps_; }
    
    /**
     * 设置载荷施加方式
     */
    void set_ramp_type(LoadRampType type) { ramp_type_ = type; }
    LoadRampType ramp_type() const { return ramp_type_; }
    
    /**
     * 设置时间范围（伪时间，用于载荷插值）
     */
    void set_time(Real start, Real end) {
        time_start_ = start;
        time_end_ = end;
    }
    Real time_start() const { return time_start_; }
    Real time_end() const { return time_end_; }
    
    // ═══ 位移边界条件 ═══
    
    /**
     * 设置位移边界条件（Dirichlet BC）
     * 
     * @param node_id 节点 ID
     * @param dof 自由度（DOF_X, DOF_Y, DOF_Z）
     * @param value_start 起始值
     * @param value_end 终止值
     */
    void set_displacement(Index node_id, int dof, Real value_start, Real value_end);
    
    /**
     * 设置位移边界条件（简化版：从 0 到 value）
     */
    void set_displacement(Index node_id, int dof, Real value) {
        set_displacement(node_id, dof, 0.0, value);
    }
    
    /**
     * 移除位移边界条件
     */
    void remove_displacement(Index node_id, int dof);
    
    /**
     * 获取位移边界条件（插值到指定时间）
     */
    std::map<std::pair<Index, int>, Real> get_displacements(Real time) const;
    
    // ═══ 力载荷 ═══
    
    /**
     * 设置节点力载荷（Neumann BC）
     * 
     * @param node_id 节点 ID
     * @param dof 自由度
     * @param value_start 起始力（N）
     * @param value_end 终止力（N）
     */
    void set_force(Index node_id, int dof, Real value_start, Real value_end);
    
    /**
     * 设置节点力载荷（简化版）
     */
    void set_force(Index node_id, int dof, Real value) {
        set_force(node_id, dof, 0.0, value);
    }
    
    /**
     * 移除力载荷
     */
    void remove_force(Index node_id, int dof);
    
    /**
     * 获取力载荷（插值到指定时间）
     */
    std::map<std::pair<Index, int>, Real> get_forces(Real time) const;
    
    // ═══ 缩放因子 ═══
    
    /**
     * 设置全局载荷缩放因子
     */
    void set_load_scale(Real scale) { load_scale_ = scale; }
    Real load_scale() const { return load_scale_; }
    
    // ═══ 工具函数 ═══
    
    /**
     * 计算给定子步的时间
     * 
     * @param substep_id 子步编号（1-based）
     * @return 对应的伪时间
     */
    Real substep_time(int substep_id) const;
    
    /**
     * 打印载荷步信息
     */
    void print() const;
    
private:
    // ═══ 数据成员 ═══
    int id_;                        // 载荷步 ID
    int num_substeps_;              // 子步数量
    LoadRampType ramp_type_;        // 载荷施加方式
    Real time_start_;               // 起始时间
    Real time_end_;                 // 终止时间
    Real load_scale_;               // 全局载荷缩放因子
    
    // 位移边界条件：(node_id, dof) -> (start_value, end_value)
    std::map<std::pair<Index, int>, std::pair<Real, Real>> displacements_;
    
    // 力载荷：(node_id, dof) -> (start_value, end_value)
    std::map<std::pair<Index, int>, std::pair<Real, Real>> forces_;
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 线性插值
     */
    Real interpolate(Real value_start, Real value_end, Real time) const;
};

/**
 * LoadStepManager - 载荷步管理器
 * 
 * 管理多个载荷步的执行
 */
class LoadStepManager {
public:
    /**
     * 添加载荷步
     */
    void add_load_step(const LoadStep& step);
    
    /**
     * 获取载荷步
     */
    const LoadStep& get_load_step(int id) const;
    LoadStep& get_load_step(int id);
    
    /**
     * 获取所有载荷步
     */
    const std::vector<LoadStep>& load_steps() const { return load_steps_; }
    
    /**
     * 载荷步数量
     */
    std::size_t num_load_steps() const { return load_steps_.size(); }
    
    /**
     * 清空所有载荷步
     */
    void clear() { load_steps_.clear(); }
    
    /**
     * 打印所有载荷步信息
     */
    void print() const;
    
private:
    std::vector<LoadStep> load_steps_;
};

}  // namespace fem
