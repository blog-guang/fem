#include "material.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace fem {
namespace constitutive {

// ═══ 默认实现：几何刚度矩阵（线性分析返回零矩阵）═══
void Material::computeGeometricStiffness(
    const Vector& /*stress*/,
    DenseMatrix& K_geo
) {
    // 默认：线性小变形分析，几何刚度为零
    K_geo.resize(strain_size_, strain_size_);
    K_geo.zero();
}

// ═══ 状态管理：默认初始化 ═══
void Material::initializeState(StateVariables& state) const {
    // 默认实现：创建与应变同尺寸的零向量
    state = StateVariables(strain_size_);
}

// ═══ 参数管理 ═══
void Material::setParameter(const std::string& name, Real value) {
    parameters_[name] = value;
}

Real Material::getParameter(const std::string& name) const {
    auto it = parameters_.find(name);
    if (it == parameters_.end()) {
        throw std::invalid_argument("Material parameter '" + name + "' not found");
    }
    return it->second;
}

std::vector<std::string> Material::getParameterNames() const {
    std::vector<std::string> names;
    names.reserve(parameters_.size());
    for (const auto& [name, _] : parameters_) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

// ═══ 工具函数 ═══
void Material::print() const {
    std::cout << "Material Type: " << typeName() << "\n";
    std::cout << "Strain Size: " << strain_size_ << "\n";
    std::cout << "Parameters:\n";
    
    auto names = getParameterNames();
    for (const auto& name : names) {
        std::cout << "  " << name << " = " << parameters_.at(name) << "\n";
    }
}

void Material::validateParameters() const {
    // 默认：无额外验证
    // 子类可重写检查参数合法性（如 E > 0, -1 < nu < 0.5 等）
}

}  // namespace constitutive
}  // namespace fem
