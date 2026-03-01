#pragma once

#include "physics/physics_base.h"
#include "physics/elasticity_unified.h"
#include "physics/heat_unified.h"
#include <string>
#include <map>
#include <memory>
#include <functional>

namespace fem {
namespace physics {

/**
 * PhysicsFactory: 物理模块工厂类
 * 
 * 提供统一的物理模块创建接口，支持：
 * - 通过字符串名称创建
 * - 参数化配置
 * - 扩展注册
 * 
 * 用法：
 * ```cpp
 * auto physics = PhysicsFactory::create("Elasticity", {
 *     {"dimension", 2}
 * }, {
 *     {"material_ptr", material}  // 需要传递材料指针
 * });
 * ```
 * 
 * 注意：对于 Elasticity，材料指针由外部管理！
 */
class PhysicsFactory {
public:
    // ═══ 参数类型 ═══
    using Parameters = std::map<std::string, Real>;
    using StringParams = std::map<std::string, std::string>;
    using PointerParams = std::map<std::string, void*>;
    
    // ═══ 核心接口 ═══
    
    /**
     * 创建物理模块实例（通用版本）
     * 
     * @param type_name 物理模块类型名称 (如 "Elasticity", "Heat")
     * @param num_params 数值参数
     * @param ptr_params 指针参数（如 material_ptr）
     * @return 物理模块指针（PhysicsBase*，调用者管理生命周期）
     * @throws std::invalid_argument 如果类型不存在或参数无效
     */
    static PhysicsBase* create(
        const std::string& type_name,
        const Parameters& num_params = {},
        const PointerParams& ptr_params = {}
    );
    
    /**
     * 创建弹性力学模块（便捷方法）
     * 
     * @param material 材料指针（外部管理）
     * @param dimension 维度 (2 or 3)
     * @return ElasticityUnified 指针
     * 
     * 示例：
     * ```cpp
     * auto* mat = new IsotropicElastic(E, nu, 2);
     * auto* physics = PhysicsFactory::createElasticity(mat, 2);
     * ```
     */
    static ElasticityUnified* createElasticity(
        constitutive::Material* material,
        int dimension
    );
    
    /**
     * 创建热传导模块（便捷方法）
     * 
     * @param conductivity 导热系数 k
     * @param source 体热源 Q
     * @return HeatConductionUnified 指针
     */
    static HeatConductionUnified* createHeat(
        Real conductivity = 1.0,
        Real source = 0.0
    );
    
    /**
     * 检查物理模块类型是否已注册
     */
    static bool isRegistered(const std::string& type_name);
    
    /**
     * 获取所有已注册的物理模块类型
     */
    static std::vector<std::string> getRegisteredTypes();
    
    /**
     * 打印所有可用物理模块类型和参数
     */
    static void printAvailableTypes();
    
    // ═══ 高级接口：自定义物理模块注册 ═══
    
    using CreatorFunc = std::function<PhysicsBase*(
        const Parameters&, const PointerParams&
    )>;
    
    /**
     * 注册新的物理模块创建器
     * 
     * @param type_name 物理模块类型名称
     * @param creator   创建函数
     */
    static void registerCreator(
        const std::string& type_name,
        CreatorFunc creator
    );

private:
    // ═══ 内部实现 ═══
    
    // 创建器注册表
    static std::map<std::string, CreatorFunc>& getRegistry();
    
    // 初始化内置物理模块类型
    static void initializeBuiltinTypes();
    
    // 辅助函数：提取参数（带默认值）
    static Real getParam(
        const Parameters& params,
        const std::string& name,
        Real default_value
    );
    
    static void* getPtrParam(
        const PointerParams& params,
        const std::string& name
    );
};

}  // namespace physics
}  // namespace fem
