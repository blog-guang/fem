#pragma once

#include "math/solver.h"
#include <string>
#include <map>
#include <memory>
#include <functional>

namespace fem {

/**
 * SolverFactory: 求解器工厂类
 * 
 * 提供统一的求解器创建接口，支持：
 * - 通过字符串名称创建
 * - 参数化配置
 * - 扩展注册
 * 
 * 用法：
 * ```cpp
 * auto solver = SolverFactory::create("PCG", {
 *     {"precond", "jacobi"},
 *     {"tol", 1e-8},
 *     {"max_iter", 1000}
 * });
 * ```
 */
class SolverFactory {
public:
    // ═══ 参数类型 ═══
    using Parameters = std::map<std::string, Real>;
    using StringParams = std::map<std::string, std::string>;
    
    // ═══ 核心接口 ═══
    
    /**
     * 创建求解器实例
     * 
     * @param type_name  求解器类型名称 (如 "CG", "PCG", "BiCGSTAB")
     * @param num_params 数值参数（tol, max_iter, omega 等）
     * @param str_params 字符串参数（precond 等）
     * @return 求解器智能指针
     * @throws std::invalid_argument 如果类型不存在
     */
    static std::unique_ptr<LinearSolver> create(
        const std::string& type_name,
        const Parameters& num_params = {},
        const StringParams& str_params = {}
    );
    
    /**
     * 创建求解器（便捷版本）
     * 
     * 示例：
     * ```cpp
     * auto solver = SolverFactory::create("CG", {{"tol", 1e-8}});
     * auto pcg = SolverFactory::create("PCG", {{"tol", 1e-10}}, {{"precond", "ilu"}});
     * ```
     */
    static std::unique_ptr<LinearSolver> create(
        const std::string& type_name,
        std::initializer_list<std::pair<const std::string, Real>> num_params,
        std::initializer_list<std::pair<const std::string, std::string>> str_params = {}
    ) {
        return create(type_name, Parameters(num_params), StringParams(str_params));
    }
    
    /**
     * 检查求解器类型是否已注册
     */
    static bool isRegistered(const std::string& type_name);
    
    /**
     * 获取所有已注册的求解器类型
     */
    static std::vector<std::string> getRegisteredTypes();
    
    /**
     * 打印所有可用求解器类型和参数
     */
    static void printAvailableTypes();
    
    // ═══ 高级接口：自定义求解器注册 ═══
    
    using CreatorFunc = std::function<std::unique_ptr<LinearSolver>(
        const Parameters&, const StringParams&
    )>;
    
    /**
     * 注册新的求解器创建器
     * 
     * @param type_name 求解器类型名称
     * @param creator   创建函数
     * 
     * 示例：
     * ```cpp
     * SolverFactory::registerCreator("MySolver",
     *     [](const auto& num_p, const auto& str_p) {
     *         auto solver = std::make_unique<MySolver>();
     *         solver->set_tol(getParam(num_p, "tol", 1e-10));
     *         return solver;
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
    
    // 初始化内置求解器类型
    static void initializeBuiltinTypes();
    
    // 辅助函数：提取参数（带默认值）
    static Real getParam(
        const Parameters& params,
        const std::string& name,
        Real default_value
    );
    
    static std::string getStrParam(
        const StringParams& params,
        const std::string& name,
        const std::string& default_value = ""
    );
};

}  // namespace fem
