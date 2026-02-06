#pragma once

#include "core/types.h"
#include <string>
#include <unordered_map>

namespace fem {

class Material {
public:
    Material(int id, const std::string& name)
        : id_(id), name_(name) {}
    
    int id() const { return id_; }
    const std::string& name() const { return name_; }
    
    // ═══ 属性管理 ═══
    void set_property(const std::string& key, Real value) {
        properties_[key] = value;
    }
    
    Real property(const std::string& key, Real default_val = 0.0) const {
        auto it = properties_.find(key);
        return (it != properties_.end()) ? it->second : default_val;
    }
    
    bool has_property(const std::string& key) const {
        return properties_.find(key) != properties_.end();
    }
    
    // ═══ 常用材料属性快捷访问 ═══
    Real E() const { return property("E"); }              // 杨氏模量
    Real nu() const { return property("nu"); }            // 泊松比
    Real rho() const { return property("rho"); }          // 密度
    Real k() const { return property("k"); }              // 热传导系数
    Real cp() const { return property("cp"); }            // 比热容
    Real alpha() const { return property("alpha"); }      // 热膨胀系数
    
    // ═══ 批量设置 ═══
    void set_elastic(Real E, Real nu) {
        set_property("E", E);
        set_property("nu", nu);
    }
    
    void set_thermal(Real k, Real rho, Real cp) {
        set_property("k", k);
        set_property("rho", rho);
        set_property("cp", cp);
    }
    
    void set_thermal_expansion(Real alpha) {
        set_property("alpha", alpha);
    }
    
    // ═══ 调试 ═══
    void print() const;
    
private:
    int id_;
    std::string name_;
    std::unordered_map<std::string, Real> properties_;
};

}  // namespace fem
