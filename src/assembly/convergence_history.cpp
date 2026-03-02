/**
 * convergence_history.cpp - Convergence History Implementation
 */

#include "assembly/convergence_history.h"
#include "core/logger.h"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace fem {

void ConvergenceHistory::add_substep(const SubstepResult& result) {
    substeps_.push_back(result);
    
    // 实时输出
    std::ostringstream oss;
    oss << "LoadStep " << result.load_step_id 
        << ", Substep " << result.substep_id
        << ", Time=" << std::fixed << std::setprecision(3) << result.time
        << ": ";
    
    if (result.converged) {
        oss << "Converged in " << result.iterations << " iterations"
            << ", Residual=" << std::scientific << std::setprecision(2) << result.residual_norm
            << ", Time=" << std::fixed << std::setprecision(3) << result.solve_time << "s";
        FEM_INFO(oss.str());
    } else {
        oss << "FAILED to converge after " << result.iterations << " iterations"
            << ", Residual=" << std::scientific << std::setprecision(2) << result.residual_norm;
        FEM_WARN(oss.str());
    }
}

int ConvergenceHistory::total_iterations() const {
    int total = 0;
    for (const auto& result : substeps_) {
        total += result.iterations;
    }
    return total;
}

Real ConvergenceHistory::total_time() const {
    Real total = 0.0;
    for (const auto& result : substeps_) {
        total += result.solve_time;
    }
    return total;
}

bool ConvergenceHistory::all_converged() const {
    for (const auto& result : substeps_) {
        if (!result.converged) {
            return false;
        }
    }
    return true;
}

void ConvergenceHistory::print_summary() const {
    FEM_INFO("========== Convergence History Summary ==========");
    FEM_INFO("Total substeps: " + std::to_string(substeps_.size()));
    FEM_INFO("Total iterations: " + std::to_string(total_iterations()));
    FEM_INFO("Total solve time: " + std::to_string(total_time()) + " s");
    FEM_INFO("All converged: " + std::string(all_converged() ? "YES" : "NO"));
    
    if (!all_converged()) {
        FEM_WARN("Failed substeps:");
        for (const auto& result : substeps_) {
            if (!result.converged) {
                FEM_WARN("  LoadStep " + std::to_string(result.load_step_id) + 
                        ", Substep " + std::to_string(result.substep_id));
            }
        }
    }
    
    FEM_INFO("=================================================");
}

void ConvergenceHistory::export_csv(const std::string& filename) const {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        FEM_ERROR("Failed to open file: " + filename);
        return;
    }
    
    // CSV header
    file << "LoadStepID,SubstepID,Time,Converged,Iterations,ResidualNorm,DisplacementNorm,SolveTime\n";
    
    // Data rows
    for (const auto& result : substeps_) {
        file << result.load_step_id << ","
             << result.substep_id << ","
             << result.time << ","
             << (result.converged ? 1 : 0) << ","
             << result.iterations << ","
             << result.residual_norm << ","
             << result.displacement_norm << ","
             << result.solve_time << "\n";
    }
    
    file.close();
    FEM_INFO("Convergence history exported to: " + filename);
}

}  // namespace fem
