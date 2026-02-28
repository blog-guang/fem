#include "math/solver_factory.h"
#include "math/cg.h"
#include "math/pcg.h"
#include "math/bicgstab.h"
#include "core/logger.h"
#include <sstream>

namespace fem {

// ═══════════════════════════════════════════════════════════
// 静态成员实现
// ═══════════════════════════════════════════════════════════

std::map<std::string, SolverFactory::CreatorFunc>& 
SolverFactory::getRegistry() {
    static std::map<std::string, CreatorFunc> registry;
    static bool initialized = false;
    
    if (!initialized) {
        initialized = true;  // 先设置，避免递归
        initializeBuiltinTypes();
    }
    
    return registry;
}

void SolverFactory::initializeBuiltinTypes() {
    FEM_INFO("SolverFactory: Initializing builtin types...");
    auto& registry = getRegistry();
    
    // ── CG（共轭梯度）──
    registry["CG"] = [](const Parameters& num_p, const StringParams& str_p) 
        -> std::unique_ptr<LinearSolver> 
    {
        auto solver = std::make_unique<CGSolver>();
        solver->set_tol(getParam(num_p, "tol", 1e-10));
        solver->set_max_iter(static_cast<std::size_t>(
            getParam(num_p, "max_iter", 1000.0)
        ));
        return solver;
    };
    
    // ── PCG（预条件共轭梯度）──
    registry["PCG"] = [](const Parameters& num_p, const StringParams& str_p) 
        -> std::unique_ptr<LinearSolver> 
    {
        std::string precond_type = getStrParam(str_p, "precond", "jacobi");
        Real omega = getParam(num_p, "omega", 1.0);
        
        auto solver = std::make_unique<PCGSolver>(precond_type, omega);
        solver->set_tol(getParam(num_p, "tol", 1e-10));
        solver->set_max_iter(static_cast<std::size_t>(
            getParam(num_p, "max_iter", 1000.0)
        ));
        return solver;
    };
    
    // ── BiCGSTAB（双共轭梯度稳定法）──
    registry["BiCGSTAB"] = [](const Parameters& num_p, const StringParams& str_p) 
        -> std::unique_ptr<LinearSolver> 
    {
        auto solver = std::make_unique<BiCGSTABSolver>();
        solver->set_tol(getParam(num_p, "tol", 1e-10));
        solver->set_max_iter(static_cast<std::size_t>(
            getParam(num_p, "max_iter", 1000.0)
        ));
        return solver;
    };
    
    FEM_INFO("SolverFactory: Registered 3 builtin solver types");
}

// ═══════════════════════════════════════════════════════════
// 核心接口
// ═══════════════════════════════════════════════════════════

std::unique_ptr<LinearSolver> SolverFactory::create(
    const std::string& type_name,
    const Parameters& num_params,
    const StringParams& str_params
) {
    auto& registry = getRegistry();
    
    auto it = registry.find(type_name);
    if (it == registry.end()) {
        std::ostringstream oss;
        oss << "Unknown solver type: '" << type_name << "'\n";
        oss << "Available types: ";
        for (const auto& [name, _] : registry) {
            oss << name << " ";
        }
        throw std::invalid_argument(oss.str());
    }
    
    try {
        auto solver = it->second(num_params, str_params);
        
        FEM_INFO("Created solver: " + type_name);
        return solver;
        
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to create solver '" << type_name << "': " << e.what();
        throw std::invalid_argument(oss.str());
    }
}

bool SolverFactory::isRegistered(const std::string& type_name) {
    auto& registry = getRegistry();
    return registry.find(type_name) != registry.end();
}

std::vector<std::string> SolverFactory::getRegisteredTypes() {
    auto& registry = getRegistry();
    std::vector<std::string> types;
    types.reserve(registry.size());
    
    for (const auto& [name, _] : registry) {
        types.push_back(name);
    }
    
    return types;
}

void SolverFactory::printAvailableTypes() {
    auto types = getRegisteredTypes();
    
    FEM_INFO("Available solver types:");
    for (const auto& name : types) {
        FEM_INFO("  - " + name);
    }
    
    FEM_INFO("");
    FEM_INFO("Parameters (all optional):");
    FEM_INFO("  Numeric:");
    FEM_INFO("    tol       : convergence tolerance (default: 1e-10)");
    FEM_INFO("    max_iter  : maximum iterations (default: 1000)");
    FEM_INFO("    omega     : relaxation parameter for SSOR (default: 1.0)");
    FEM_INFO("  String:");
    FEM_INFO("    precond   : preconditioner type (PCG only)");
    FEM_INFO("                jacobi, ssor, ilu, amg (default: jacobi)");
}

void SolverFactory::registerCreator(
    const std::string& type_name,
    CreatorFunc creator
) {
    auto& registry = getRegistry();
    
    if (registry.find(type_name) != registry.end()) {
        FEM_WARN("Overwriting existing solver type: " + type_name);
    }
    
    registry[type_name] = creator;
    FEM_INFO("Registered solver type: " + type_name);
}

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

Real SolverFactory::getParam(
    const Parameters& params,
    const std::string& name,
    Real default_value
) {
    auto it = params.find(name);
    return (it != params.end()) ? it->second : default_value;
}

std::string SolverFactory::getStrParam(
    const StringParams& params,
    const std::string& name,
    const std::string& default_value
) {
    auto it = params.find(name);
    return (it != params.end()) ? it->second : default_value;
}

}  // namespace fem
