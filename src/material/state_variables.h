#pragma once

#include "core/types.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include <map>
#include <string>
#include <iosfwd>

namespace fem {
namespace constitutive {

/**
 * StateVariables: 材料内部状态变量容器
 * 
 * 存储塑性应变、损伤、硬化变量等历史相关量。
 * 支持序列化/反序列化用于检查点重启。
 */
class StateVariables {
public:
    // ═══ 构造 ═══
    StateVariables() = default;
    explicit StateVariables(std::size_t strain_size);
    
    // ═══ 常用状态变量 (预定义成员) ═══
    Vector plastic_strain;           // 塑性应变张量 (Voigt记号)
    Real   equiv_plastic_strain{0.0}; // 等效塑性应变 (累积)
    Real   damage{0.0};              // 损伤变量 [0,1]
    Vector back_stress;              // 背应力张量 (运动硬化)
    
    // ═══ 通用扩展变量 (map存储) ═══
    std::map<std::string, Real>   scalar_vars;  // 标量变量
    std::map<std::string, Vector> tensor_vars;  // 张量变量
    
    // ═══ 访问接口 ═══
    // 设置/获取标量变量
    void setScalar(const std::string& name, Real value);
    Real getScalar(const std::string& name, Real default_val = 0.0) const;
    
    // 设置/获取张量变量
    void setTensor(const std::string& name, const Vector& value);
    Vector getTensor(const std::string& name) const;
    bool hasTensor(const std::string& name) const;
    
    // ═══ 工具函数 ═══
    // 重置所有状态为初始值
    void reset();
    
    // 深拷贝
    StateVariables clone() const;
    
    // 序列化/反序列化 (用于检查点)
    void serialize(std::ostream& os) const;
    void deserialize(std::istream& is);
    
    // 调试输出
    void print(const std::string& prefix = "") const;
};

}  // namespace constitutive
}  // namespace fem
