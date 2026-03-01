#include "physics/physics_factory.h"
#include "core/logger.h"
#include <sstream>

namespace fem {
namespace physics {

// ═══════════════════════════════════════════════════════════
// 静态成员实现
// ═══════════════════════════════════════════════════════════

std::map<std::string, PhysicsFactory::CreatorFunc>& 
PhysicsFactory::getRegistry() {
    static std::map<std::string, CreatorFunc> registry;
    static bool initialized = false;
    
    if (!initialized) {
        initialized = true;  // 先设置，避免递归
        initializeBuiltinTypes();
    }
    
    return registry;
}

void PhysicsFactory::initializeBuiltinTypes() {
    FEM_INFO("PhysicsFactory: Initializing builtin types...");
    auto& registry = getRegistry();
    
    // ── Elasticity ──
    registry["Elasticity"] = [](const Parameters& num_p, const PointerParams& ptr_p) 
        -> PhysicsBase* 
    {
        auto* material = static_cast<constitutive::Material*>(
            getPtrParam(ptr_p, "material_ptr")
        );
        int dimension = static_cast<int>(getParam(num_p, "dimension", 3.0));
        
        return new ElasticityUnified(material, dimension);
    };
    
    // ── Heat ──
    registry["Heat"] = [](const Parameters& num_p, const PointerParams& ptr_p) 
        -> PhysicsBase* 
    {
        Real k = getParam(num_p, "conductivity", 1.0);
        Real Q = getParam(num_p, "source", 0.0);
        
        return new HeatConductionUnified(k, Q);
    };
    
    FEM_INFO("PhysicsFactory: Registered 2 builtin physics types");
}

// ═══════════════════════════════════════════════════════════
// 核心接口
// ═══════════════════════════════════════════════════════════

PhysicsBase* PhysicsFactory::create(
    const std::string& type_name,
    const Parameters& num_params,
    const PointerParams& ptr_params
) {
    auto& registry = getRegistry();
    
    auto it = registry.find(type_name);
    if (it == registry.end()) {
        std::ostringstream oss;
        oss << "Unknown physics type: '" << type_name << "'\n";
        oss << "Available types: ";
        for (const auto& [name, _] : registry) {
            oss << name << " ";
        }
        throw std::invalid_argument(oss.str());
    }
    
    try {
        PhysicsBase* physics = it->second(num_params, ptr_params);
        
        FEM_INFO("Created physics module: " + type_name);
        return physics;
        
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to create physics module '" << type_name << "': " << e.what();
        throw std::invalid_argument(oss.str());
    }
}

ElasticityUnified* PhysicsFactory::createElasticity(
    constitutive::Material* material,
    int dimension
) {
    return new ElasticityUnified(material, dimension);
}

HeatConductionUnified* PhysicsFactory::createHeat(
    Real conductivity,
    Real source
) {
    return new HeatConductionUnified(conductivity, source);
}

bool PhysicsFactory::isRegistered(const std::string& type_name) {
    auto& registry = getRegistry();
    return registry.find(type_name) != registry.end();
}

std::vector<std::string> PhysicsFactory::getRegisteredTypes() {
    auto& registry = getRegistry();
    std::vector<std::string> types;
    types.reserve(registry.size());
    
    for (const auto& [name, _] : registry) {
        types.push_back(name);
    }
    
    return types;
}

void PhysicsFactory::printAvailableTypes() {
    auto types = getRegisteredTypes();
    
    FEM_INFO("Available physics types:");
    for (const auto& name : types) {
        FEM_INFO("  - " + name);
    }
    
    FEM_INFO("");
    FEM_INFO("Parameters:");
    FEM_INFO("  Elasticity:");
    FEM_INFO("    material_ptr : Material* (pointer, required)");
    FEM_INFO("    dimension    : 2 or 3 (default: 3)");
    FEM_INFO("  Heat:");
    FEM_INFO("    conductivity : Real (default: 1.0)");
    FEM_INFO("    source       : Real (default: 0.0)");
}

void PhysicsFactory::registerCreator(
    const std::string& type_name,
    CreatorFunc creator
) {
    auto& registry = getRegistry();
    
    if (registry.find(type_name) != registry.end()) {
        FEM_WARN("Overwriting existing physics type: " + type_name);
    }
    
    registry[type_name] = creator;
    FEM_INFO("Registered physics type: " + type_name);
}

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

Real PhysicsFactory::getParam(
    const Parameters& params,
    const std::string& name,
    Real default_value
) {
    auto it = params.find(name);
    return (it != params.end()) ? it->second : default_value;
}

void* PhysicsFactory::getPtrParam(
    const PointerParams& params,
    const std::string& name
) {
    auto it = params.find(name);
    if (it == params.end()) {
        throw std::invalid_argument(
            "Missing required pointer parameter: '" + name + "'"
        );
    }
    return it->second;
}

}  // namespace physics
}  // namespace fem
