/**
 * data_manager.h - 场数据管理器
 * 
 * 管理网格中的所有场数据
 */

#pragma once

#include "data/base_data.h"
#include "data/field_data.h"
#include "core/logger.h"
#include <map>
#include <memory>
#include <vector>
#include <iostream>

namespace fem {
namespace data {

/**
 * DataManager: 场数据管理器
 * 
 * 功能：
 * - 注册和管理多个场数据
 * - 按名称查找场
 * - 批量操作（清空、重置、调整大小）
 * - 数据导出和导入
 * 
 * 使用示例：
 * ```cpp
 * DataManager manager;
 * 
 * // 注册场
 * manager.add_field<RealData>("temperature", DataLocation::Node, n_nodes);
 * manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes);
 * manager.add_field<TensorData>("stress", DataLocation::GaussPoint, n_elems * n_gp);
 * 
 * // 访问场
 * auto* temp = manager.get_field<RealData>("temperature");
 * temp->set(node_id, 300.0);
 * Real T = temp->get(node_id);
 * 
 * // 批量操作
 * manager.resize_all(DataLocation::Node, new_n_nodes);
 * manager.clear_all();
 * ```
 */
class DataManager {
public:
    DataManager() = default;
    ~DataManager() = default;

    // ═══ 场注册 ═══
    
    /**
     * 添加场数据
     * @tparam T 场数据类型（RealData, VectorData等）
     * @param name 场名称
     * @param location 数据位置
     * @param size 初始大小
     * @return 场数据指针
     */
    template<typename T>
    T* add_field(const std::string& name, DataLocation location, Index size = 0) {
        if (has_field(name)) {
            FEM_WARN("DataManager: field '" + name + "' already exists");
            return get_field<T>(name);
        }
        
        auto field = std::make_unique<T>(name, location, size);
        T* ptr = field.get();
        fields_[name] = std::move(field);
        return ptr;
    }
    
    /**
     * 添加场数据（带默认值）
     */
    template<typename T, typename ValueType>
    T* add_field(const std::string& name, DataLocation location, 
                 Index size, const ValueType& default_value) {
        if (has_field(name)) {
            FEM_WARN("DataManager: field '" + name + "' already exists");
            return get_field<T>(name);
        }
        
        auto field = std::make_unique<T>(name, location, size, default_value);
        T* ptr = field.get();
        fields_[name] = std::move(field);
        return ptr;
    }
    
    /**
     * 移除场数据
     */
    bool remove_field(const std::string& name) {
        return fields_.erase(name) > 0;
    }
    
    /**
     * 检查场是否存在
     */
    bool has_field(const std::string& name) const {
        return fields_.find(name) != fields_.end();
    }

    // ═══ 场访问 ═══
    
    /**
     * 获取场数据（类型安全）
     * @tparam T 场数据类型
     * @param name 场名称
     * @return 场数据指针（失败返回 nullptr）
     */
    template<typename T>
    T* get_field(const std::string& name) {
        auto it = fields_.find(name);
        if (it == fields_.end()) {
            FEM_WARN("DataManager: field '" + name + "' not found");
            return nullptr;
        }
        return dynamic_cast<T*>(it->second.get());
    }
    
    /**
     * 获取场数据（常量版本）
     */
    template<typename T>
    const T* get_field(const std::string& name) const {
        auto it = fields_.find(name);
        if (it == fields_.end()) {
            return nullptr;
        }
        return dynamic_cast<const T*>(it->second.get());
    }
    
    /**
     * 获取基类指针（用于通用操作）
     */
    BaseData* get_field_base(const std::string& name) {
        auto it = fields_.find(name);
        if (it == fields_.end()) {
            return nullptr;
        }
        return it->second.get();
    }
    
    /**
     * 获取所有场名称
     */
    std::vector<std::string> get_field_names() const {
        std::vector<std::string> names;
        names.reserve(fields_.size());
        for (const auto& pair : fields_) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    /**
     * 获取指定位置的所有场名称
     */
    std::vector<std::string> get_field_names(DataLocation location) const {
        std::vector<std::string> names;
        for (const auto& pair : fields_) {
            if (pair.second->location() == location) {
                names.push_back(pair.first);
            }
        }
        return names;
    }

    // ═══ 批量操作 ═══
    
    /**
     * 调整指定位置所有场的大小
     */
    void resize_all(DataLocation location, Index new_size) {
        for (auto& pair : fields_) {
            if (pair.second->location() == location) {
                pair.second->resize(new_size);
            }
        }
    }
    
    /**
     * 清空所有场
     */
    void clear_all() {
        for (auto& pair : fields_) {
            pair.second->clear();
        }
    }
    
    /**
     * 重置所有场为默认值
     */
    void reset_all() {
        for (auto& pair : fields_) {
            pair.second->reset();
        }
    }
    
    /**
     * 重置指定位置的所有场
     */
    void reset_all(DataLocation location) {
        for (auto& pair : fields_) {
            if (pair.second->location() == location) {
                pair.second->reset();
            }
        }
    }

    // ═══ 统计信息 ═══
    
    /**
     * 获取场数量
     */
    size_t num_fields() const { return fields_.size(); }
    
    /**
     * 获取指定位置的场数量
     */
    size_t num_fields(DataLocation location) const {
        size_t count = 0;
        for (const auto& pair : fields_) {
            if (pair.second->location() == location) {
                ++count;
            }
        }
        return count;
    }
    
    /**
     * 打印所有场的信息
     */
    void print_info() const {
        std::cout << "DataManager: " << fields_.size() << " fields\n";
        std::cout << "───────────────────────────────────────\n";
        for (const auto& pair : fields_) {
            pair.second->print_info();
            std::cout << "───────────────────────────────────────\n";
        }
    }

private:
    std::map<std::string, std::unique_ptr<BaseData>> fields_;
};

} // namespace data
} // namespace fem
