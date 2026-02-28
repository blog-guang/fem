#pragma once

#include "material/material.h"
#include <string>
#include <map>
#include <functional>
#include <stdexcept>

namespace fem {
namespace constitutive {

/**
 * MaterialFactory: 材料工厂类
 * 
 * 提供统一的材料创建接口，支持：
 * - 通过字符串名称创建材料
 * - 参数化创建
 * - 材料类型注册（可扩展）
 * 
 * 用法：
 * ```cpp
 * auto mat = MaterialFactory::create("IsotropicElastic", {
 *     {"E", 200e9},
 *     {"nu", 0.3},
 *     {"dimension", 3}
 * });
 * ```
 */
class MaterialFactory {
public:
    // ═══ 参数类型 ═══
    using Parameters = std::map<std::string, Real>;
    
    // ═══ 核心接口 ═══
    
    /**
     * 创建材料实例
     * 
     * @param type_name 材料类型名称 (如 "IsotropicElastic", "J2Plasticity")
     * @param params    材料参数
     * @return 材料智能指针
     * @throws std::invalid_argument 如果类型不存在或参数无效
     */
    static MaterialPtr create(
        const std::string& type_name,
        const Parameters& params
    );
    
    /**
     * 创建材料实例（便捷版本，支持初始化列表）
     * 
     * @param type_name 材料类型名称
     * @param params    参数初始化列表
     * 
     * 示例：
     * ```cpp
     * auto mat = MaterialFactory::create("IsotropicElastic", {
     *     {"E", 200e9}, {"nu", 0.3}, {"dimension", 3}
     * });
     * ```
     */
    static MaterialPtr create(
        const std::string& type_name,
        std::initializer_list<std::pair<const std::string, Real>> params
    ) {
        return create(type_name, Parameters(params));
    }
    
    /**
     * 检查材料类型是否已注册
     */
    static bool isRegistered(const std::string& type_name);
    
    /**
     * 获取所有已注册的材料类型名称
     */
    static std::vector<std::string> getRegisteredTypes();
    
    /**
     * 打印所有可用材料类型和所需参数
     */
    static void printAvailableTypes();
    
    // ═══ 高级接口：自定义材料注册 ═══
    
    using CreatorFunc = std::function<MaterialPtr(const Parameters&)>;
    
    /**
     * 注册新的材料创建器
     * 
     * @param type_name 材料类型名称
     * @param creator   创建函数
     * 
     * 示例：
     * ```cpp
     * MaterialFactory::registerCreator("MyMaterial", 
     *     [](const auto& p) -> MaterialPtr {
     *         return std::make_shared<MyMaterial>(
     *             p.at("param1"), p.at("param2")
     *         );
     *     }
     * );
     * ```
     */
    static void registerCreator(
        const std::string& type_name,
        CreatorFunc creator
    );

private:
    // ═══ 内部实现 ═══
    
    // 创建器注册表
    static std::map<std::string, CreatorFunc>& getRegistry();
    
    // 初始化内置材料类型
    static void initializeBuiltinTypes();
    
    // 辅助函数：提取参数（带默认值）
    static Real getParam(
        const Parameters& params,
        const std::string& name,
        Real default_value = 0.0
    );
    
    // 辅助函数：提取必需参数（不存在则抛异常）
    static Real requireParam(
        const Parameters& params,
        const std::string& name
    );
};

}  // namespace constitutive
}  // namespace fem
