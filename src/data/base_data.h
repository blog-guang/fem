/**
 * base_data.h - 场数据基类
 * 
 * 定义场数据的通用接口
 */

#pragma once

#include "core/types.h"
#include "data/data_location.h"
#include <string>
#include <memory>

namespace fem {
namespace data {

/**
 * BaseData: 场数据基类
 * 
 * 抽象基类，定义场数据的通用接口：
 * - 数据存储位置（Node, Element, Face, GaussPoint等）
 * - 数据大小管理
 * - 清空和重置
 * - 克隆
 * 
 * 子类实现具体的数据类型（Real, Int, Vector, Tensor等）
 */
class BaseData {
public:
    /**
     * 构造函数
     * @param name 场名称
     * @param location 数据存储位置
     */
    BaseData(const std::string& name, DataLocation location)
        : name_(name), location_(location), size_(0) {}
    
    virtual ~BaseData() = default;

    // ═══ 基本信息访问 ═══
    
    /**
     * 获取场名称
     */
    const std::string& name() const { return name_; }
    
    /**
     * 获取数据存储位置
     */
    DataLocation location() const { return location_; }
    
    /**
     * 获取数据大小（存储单元数量）
     */
    Index size() const { return size_; }
    
    /**
     * 是否为空
     */
    bool empty() const { return size_ == 0; }
    
    /**
     * 获取数据类型名称（用于调试和序列化）
     */
    virtual std::string type_name() const = 0;

    // ═══ 数据管理 ═══
    
    /**
     * 调整数据大小
     * @param new_size 新的大小
     */
    virtual void resize(Index new_size) = 0;
    
    /**
     * 清空数据（释放内存）
     */
    virtual void clear() = 0;
    
    /**
     * 重置数据为默认值（保持大小）
     */
    virtual void reset() = 0;
    
    /**
     * 克隆数据（深拷贝）
     */
    virtual std::unique_ptr<BaseData> clone() const = 0;

    // ═══ 序列化支持（可选）═══
    
    /**
     * 打印信息（用于调试）
     */
    virtual void print_info() const;

protected:
    std::string name_;        ///< 场名称
    DataLocation location_;   ///< 数据存储位置
    Index size_;              ///< 数据大小
};

} // namespace data
} // namespace fem
