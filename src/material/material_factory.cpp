#include "material/material_factory.h"
#include "material/isotropic_elastic.h"
#include "material/orthotropic_elastic.h"
#include "material/j2_plasticity.h"
#include "material/j2_plasticity_kinematic.h"
#include "material/neo_hookean.h"
#include "core/logger.h"
#include <sstream>

namespace fem {
namespace constitutive {

// ═══════════════════════════════════════════════════════════
// 静态成员实现
// ═══════════════════════════════════════════════════════════

std::map<std::string, MaterialFactory::CreatorFunc>& 
MaterialFactory::getRegistry() {
    static std::map<std::string, CreatorFunc> registry;
    static bool initialized = false;
    
    if (!initialized) {
        initialized = true;  // 先设置标志，避免递归
        initializeBuiltinTypes();
    }
    
    return registry;
}

void MaterialFactory::initializeBuiltinTypes() {
    FEM_INFO("MaterialFactory: Initializing builtin types...");
    auto& registry = getRegistry();
    FEM_INFO("MaterialFactory: Registry obtained");
    
    // ── IsotropicElastic ──
    registry["IsotropicElastic"] = [](const Parameters& p) -> MaterialPtr {
        Real E = requireParam(p, "E");
        Real nu = requireParam(p, "nu");
        int dim = static_cast<int>(getParam(p, "dimension", 3.0));
        bool plane_stress = getParam(p, "plane_stress", 1.0) > 0.5;
        
        return std::make_shared<IsotropicElastic>(E, nu, dim, plane_stress);
    };
    
    // ── OrthotropicElastic ──
    registry["OrthotropicElastic"] = [](const Parameters& p) -> MaterialPtr {
        Real E1 = requireParam(p, "E1");
        Real E2 = requireParam(p, "E2");
        Real nu12 = requireParam(p, "nu12");
        Real G12 = requireParam(p, "G12");
        
        // 可选参数（3D）
        Real E3 = getParam(p, "E3", E1);
        Real nu13 = getParam(p, "nu13", nu12);
        Real nu23 = getParam(p, "nu23", nu12);
        Real G13 = getParam(p, "G13", G12);
        Real G23 = getParam(p, "G23", G12);
        
        int dim = static_cast<int>(getParam(p, "dimension", 2.0));
        
        if (dim == 2) {
            return std::make_shared<OrthotropicElastic>(
                E1, E2, nu12, G12
            );
        } else {
            return std::make_shared<OrthotropicElastic>(
                E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
            );
        }
    };
    
    // ── J2Plasticity（等向硬化）──
    registry["J2Plasticity"] = [](const Parameters& p) -> MaterialPtr {
        Real E = requireParam(p, "E");
        Real nu = requireParam(p, "nu");
        Real sigma_y = requireParam(p, "sigma_y");
        Real H = getParam(p, "H", 0.0);  // 默认理想塑性
        int dim = static_cast<int>(getParam(p, "dimension", 3.0));
        
        return std::make_shared<J2Plasticity>(E, nu, sigma_y, H, dim);
    };
    
    // ── J2PlasticityKinematic（随动硬化）──
    registry["J2PlasticityKinematic"] = [](const Parameters& p) -> MaterialPtr {
        Real E = requireParam(p, "E");
        Real nu = requireParam(p, "nu");
        Real sigma_y = requireParam(p, "sigma_y");
        Real H_iso = getParam(p, "H_iso", 0.0);   // 等向硬化
        Real H_kin = getParam(p, "H_kin", 0.0);   // 随动硬化
        int dim = static_cast<int>(getParam(p, "dimension", 3.0));
        
        return std::make_shared<J2PlasticityKinematic>(
            E, nu, sigma_y, H_iso, H_kin, dim
        );
    };
    
    // ── NeoHookean（超弹性）──
    registry["NeoHookean"] = [](const Parameters& p) -> MaterialPtr {
        int dim = static_cast<int>(getParam(p, "dimension", 3.0));
        
        // 支持两种参数输入方式
        if (p.find("C10") != p.end() && p.find("D1") != p.end()) {
            // 直接使用 C10, D1 参数
            Real C10 = requireParam(p, "C10");
            Real D1 = requireParam(p, "D1");
            return std::make_shared<NeoHookean>(C10, D1, dim);
        } else {
            // 从工程参数 E, nu 转换
            Real E = requireParam(p, "E");
            Real nu = requireParam(p, "nu");
            return std::make_shared<NeoHookean>(E, nu, dim, true);
        }
    };
}

// ═══════════════════════════════════════════════════════════
// 核心接口
// ═══════════════════════════════════════════════════════════

MaterialPtr MaterialFactory::create(
    const std::string& type_name,
    const Parameters& params
) {
    auto& registry = getRegistry();
    
    auto it = registry.find(type_name);
    if (it == registry.end()) {
        std::ostringstream oss;
        oss << "Unknown material type: '" << type_name << "'\n";
        oss << "Available types: ";
        for (const auto& [name, _] : registry) {
            oss << name << " ";
        }
        throw std::invalid_argument(oss.str());
    }
    
    try {
        MaterialPtr mat = it->second(params);
        
        // 验证参数合法性
        mat->validateParameters();
        
        FEM_INFO("Created material: " + type_name);
        return mat;
        
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to create material '" << type_name << "': " << e.what();
        throw std::invalid_argument(oss.str());
    }
}

bool MaterialFactory::isRegistered(const std::string& type_name) {
    auto& registry = getRegistry();
    return registry.find(type_name) != registry.end();
}

std::vector<std::string> MaterialFactory::getRegisteredTypes() {
    auto& registry = getRegistry();
    std::vector<std::string> types;
    types.reserve(registry.size());
    
    for (const auto& [name, _] : registry) {
        types.push_back(name);
    }
    
    return types;
}

void MaterialFactory::printAvailableTypes() {
    auto types = getRegisteredTypes();
    
    FEM_INFO("Available material types:");
    for (const auto& name : types) {
        FEM_INFO("  - " + name);
    }
    
    FEM_INFO("");
    FEM_INFO("Required parameters:");
    FEM_INFO("  IsotropicElastic: E, nu, [dimension=3], [plane_stress=1]");
    FEM_INFO("  OrthotropicElastic: E1, E2, nu12, G12, [dimension=2]");
    FEM_INFO("  J2Plasticity: E, nu, sigma_y, [H=0], [dimension=3]");
    FEM_INFO("  J2PlasticityKinematic: E, nu, sigma_y, [H_iso=0], [H_kin=0], [dimension=3]");
}

void MaterialFactory::registerCreator(
    const std::string& type_name,
    CreatorFunc creator
) {
    auto& registry = getRegistry();
    
    if (registry.find(type_name) != registry.end()) {
        FEM_WARN("Overwriting existing material type: " + type_name);
    }
    
    registry[type_name] = creator;
    FEM_INFO("Registered material type: " + type_name);
}

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

Real MaterialFactory::getParam(
    const Parameters& params,
    const std::string& name,
    Real default_value
) {
    auto it = params.find(name);
    return (it != params.end()) ? it->second : default_value;
}

Real MaterialFactory::requireParam(
    const Parameters& params,
    const std::string& name
) {
    auto it = params.find(name);
    if (it == params.end()) {
        throw std::invalid_argument(
            "Missing required parameter: '" + name + "'"
        );
    }
    return it->second;
}

}  // namespace constitutive
}  // namespace fem
