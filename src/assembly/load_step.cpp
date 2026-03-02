/**
 * load_step.cpp - Load Step Implementation
 */

#include "assembly/load_step.h"
#include "core/logger.h"
#include <sstream>
#include <stdexcept>

namespace fem {

// ═══════════════════════════════════════════════════════════
// LoadStep 实现
// ═══════════════════════════════════════════════════════════

void LoadStep::set_displacement(Index node_id, int dof, Real value_start, Real value_end) {
    auto key = std::make_pair(node_id, dof);
    displacements_[key] = std::make_pair(value_start, value_end);
}

void LoadStep::remove_displacement(Index node_id, int dof) {
    auto key = std::make_pair(node_id, dof);
    displacements_.erase(key);
}

std::map<std::pair<Index, int>, Real> LoadStep::get_displacements(Real time) const {
    std::map<std::pair<Index, int>, Real> result;
    
    for (const auto& [key, values] : displacements_) {
        Real value = interpolate(values.first, values.second, time);
        result[key] = value * load_scale_;
    }
    
    return result;
}

void LoadStep::set_force(Index node_id, int dof, Real value_start, Real value_end) {
    auto key = std::make_pair(node_id, dof);
    forces_[key] = std::make_pair(value_start, value_end);
}

void LoadStep::remove_force(Index node_id, int dof) {
    auto key = std::make_pair(node_id, dof);
    forces_.erase(key);
}

std::map<std::pair<Index, int>, Real> LoadStep::get_forces(Real time) const {
    std::map<std::pair<Index, int>, Real> result;
    
    for (const auto& [key, values] : forces_) {
        Real value = interpolate(values.first, values.second, time);
        result[key] = value * load_scale_;
    }
    
    return result;
}

Real LoadStep::substep_time(int substep_id) const {
    if (substep_id < 1 || substep_id > num_substeps_) {
        throw std::out_of_range("Substep ID out of range");
    }
    
    // 线性插值时间
    Real lambda = Real(substep_id) / Real(num_substeps_);
    return time_start_ + lambda * (time_end_ - time_start_);
}

Real LoadStep::interpolate(Real value_start, Real value_end, Real time) const {
    // 归一化时间 [time_start, time_end] -> [0, 1]
    Real time_range = time_end_ - time_start_;
    Real lambda = 0.0;
    
    if (time_range > 1e-15) {
        lambda = (time - time_start_) / time_range;
        lambda = std::max(0.0, std::min(1.0, lambda));  // 限制在 [0,1]
    }
    
    // 根据载荷施加方式插值
    if (ramp_type_ == LoadRampType::RAMPED) {
        // 线性递增
        return value_start + lambda * (value_end - value_start);
    } else {
        // 阶跃：lambda > 0 时直接用终止值
        return (lambda > 0.0) ? value_end : value_start;
    }
}

void LoadStep::print() const {
    std::ostringstream oss;
    oss << "LoadStep " << id_ << ":\n";
    oss << "  Substeps: " << num_substeps_ << "\n";
    oss << "  Ramp: " << (ramp_type_ == LoadRampType::RAMPED ? "Ramped" : "Stepped") << "\n";
    oss << "  Time: [" << time_start_ << ", " << time_end_ << "]\n";
    oss << "  Load Scale: " << load_scale_ << "\n";
    oss << "  Displacements: " << displacements_.size() << " BCs\n";
    oss << "  Forces: " << forces_.size() << " loads\n";
    
    FEM_INFO(oss.str());
}

// ═══════════════════════════════════════════════════════════
// LoadStepManager 实现
// ═══════════════════════════════════════════════════════════

void LoadStepManager::add_load_step(const LoadStep& step) {
    load_steps_.push_back(step);
}

const LoadStep& LoadStepManager::get_load_step(int id) const {
    for (const auto& step : load_steps_) {
        if (step.id() == id) {
            return step;
        }
    }
    throw std::out_of_range("LoadStep ID not found: " + std::to_string(id));
}

LoadStep& LoadStepManager::get_load_step(int id) {
    for (auto& step : load_steps_) {
        if (step.id() == id) {
            return step;
        }
    }
    throw std::out_of_range("LoadStep ID not found: " + std::to_string(id));
}

void LoadStepManager::print() const {
    FEM_INFO("LoadStepManager: " + std::to_string(load_steps_.size()) + " load steps");
    for (const auto& step : load_steps_) {
        step.print();
    }
}

}  // namespace fem
