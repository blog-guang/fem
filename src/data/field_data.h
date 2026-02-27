/**
 * field_data.h - 模板化场数据类
 * 
 * 提供类型安全的场数据存储和访问
 */

#pragma once

#include "data/base_data.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "core/logger.h"
#include <vector>
#include <stdexcept>
#include <sstream>
#include <type_traits>

namespace fem {
namespace data {

/**
 * FieldData<T>: 类型化场数据
 * 
 * 模板参数 T 可以是：
 * - Real: 标量场（温度、压力、密度等）
 * - int: 整型场（材料ID、标签等）
 * - bool: 布尔场（激活标志等）
 * - std::string: 字符串场（名称等）
 * - Vector: 矢量场（位移、速度等）
 * - DenseMatrix: 张量场（应力、应变等）
 * 
 * 数据存储方式：
 * - 对于 Node 位置：data_[node_id] = value
 * - 对于 Element 位置：data_[elem_id] = value
 * - 对于 GaussPoint 位置：data_[elem_id * n_gp + gp_idx] = value
 */
template<typename T>
class FieldData : public BaseData {
public:
    /**
     * 构造函数
     * @param name 场名称
     * @param location 数据存储位置
     * @param initial_size 初始大小（默认为0）
     */
    FieldData(const std::string& name, 
              DataLocation location,
              Index initial_size = 0)
        : BaseData(name, location) {
        if (initial_size > 0) {
            resize(initial_size);
        }
    }
    
    /**
     * 构造函数（带默认值）
     * @param name 场名称
     * @param location 数据存储位置
     * @param initial_size 初始大小
     * @param default_value 默认值
     */
    FieldData(const std::string& name,
              DataLocation location,
              Index initial_size,
              const T& default_value)
        : BaseData(name, location), default_value_(default_value) {
        resize(initial_size);
        reset();
    }

    // ═══ 数据访问 ═══
    
    /**
     * 设置数据（单个索引）
     * @param index 索引（node_id, elem_id, 或 gauss_point_id）
     * @param value 值
     */
    void set(Index index, const T& value) {
        check_index(index);
        data_[index] = value;
    }
    
    /**
     * 获取数据（单个索引）
     * @param index 索引
     * @return 值的常引用
     */
    const T& get(Index index) const {
        check_index(index);
        return data_[index];
    }
    
    /**
     * 获取数据（单个索引，可修改）
     * @param index 索引
     * @return 值的引用
     */
    T& get(Index index) {
        check_index(index);
        return data_[index];
    }
    
    /**
     * 设置高斯点数据
     * @param elem_id 单元ID
     * @param gp_index 高斯点索引（0-based）
     * @param n_gp_per_elem 每单元高斯点数量
     * @param value 值
     */
    void set_gauss_point(Index elem_id, Index gp_index, Index n_gp_per_elem, const T& value) {
        if (location_ != DataLocation::GaussPoint) {
            FEM_WARN("FieldData::set_gauss_point called on non-GaussPoint field");
        }
        Index index = elem_id * n_gp_per_elem + gp_index;
        set(index, value);
    }
    
    /**
     * 获取高斯点数据
     * @param elem_id 单元ID
     * @param gp_index 高斯点索引
     * @param n_gp_per_elem 每单元高斯点数量
     * @return 值的常引用
     */
    const T& get_gauss_point(Index elem_id, Index gp_index, Index n_gp_per_elem) const {
        Index index = elem_id * n_gp_per_elem + gp_index;
        return get(index);
    }
    
    /**
     * 批量设置数据
     * @param values 值数组
     */
    void set_all(const std::vector<T>& values) {
        if (values.size() != size_) {
            FEM_WARN("FieldData::set_all size mismatch");
            resize(values.size());
        }
        data_ = values;
    }
    
    /**
     * 获取所有数据
     */
    const std::vector<T>& data() const { return data_; }
    
    /**
     * 获取所有数据（可修改）
     */
    std::vector<T>& data() { return data_; }
    
    /**
     * 运算符重载：直接索引
     */
    const T& operator[](Index index) const { return get(index); }
    T& operator[](Index index) { return get(index); }

    // ═══ BaseData 接口实现 ═══
    
    void resize(Index new_size) override {
        data_.resize(new_size, default_value_);
        size_ = new_size;
    }
    
    void clear() override {
        data_.clear();
        size_ = 0;
    }
    
    void reset() override {
        std::fill(data_.begin(), data_.end(), default_value_);
    }
    
    std::string type_name() const override {
        return get_type_name();
    }
    
    std::unique_ptr<BaseData> clone() const override {
        auto cloned = std::make_unique<FieldData<T>>(name_, location_, size_);
        cloned->data_ = data_;
        cloned->default_value_ = default_value_;
        return cloned;
    }

    // ═══ 统计信息 ═══
    
    /**
     * 获取默认值
     */
    const T& default_value() const { return default_value_; }
    
    /**
     * 设置默认值
     */
    void set_default_value(const T& value) { default_value_ = value; }

private:
    std::vector<T> data_;       ///< 数据存储
    T default_value_{};         ///< 默认值
    
    /**
     * 检查索引有效性
     */
    void check_index(Index index) const {
        if (index >= size_) {
            std::ostringstream oss;
            oss << "FieldData '" << name_ << "' index out of range: "
                << index << " >= " << size_;
            throw std::out_of_range(oss.str());
        }
    }
    
    /**
     * 获取类型名称（用于调试）
     */
    static std::string get_type_name() {
        if (std::is_same<T, Real>::value) return "Real";
        if (std::is_same<T, int>::value) return "Int";
        if (std::is_same<T, bool>::value) return "Bool";
        if (std::is_same<T, std::string>::value) return "String";
        if (std::is_same<T, Vector>::value) return "Vector";
        if (std::is_same<T, DenseMatrix>::value) return "Tensor";
        return "Unknown";
    }
};

// ═══════════════════════════════════════════════════════════
// 类型别名（方便使用）
// ═══════════════════════════════════════════════════════════

// ═══════════════════════════════════════════════════════════
// BoolData 特化（处理 std::vector<bool> 的特殊性）
// ═══════════════════════════════════════════════════════════

/**
 * BoolData 特化版本
 * 
 * 由于 std::vector<bool> 的特殊实现，需要特化处理
 */
template<>
class FieldData<bool> : public BaseData {
public:
    FieldData(const std::string& name, DataLocation location, Index initial_size = 0)
        : BaseData(name, location), default_value_(false) {
        if (initial_size > 0) {
            resize(initial_size);
        }
    }
    
    FieldData(const std::string& name, DataLocation location, 
              Index initial_size, bool default_value)
        : BaseData(name, location), default_value_(default_value) {
        resize(initial_size);
        reset();
    }
    
    // 设置值
    void set(Index index, bool value) {
        check_index(index);
        data_[index] = value;
    }
    
    // 获取值（返回值而非引用）
    bool get(Index index) const {
        check_index(index);
        return data_[index];
    }
    
    // 运算符重载（返回值）
    bool operator[](Index index) const { return get(index); }
    
    // 批量设置
    void set_all(const std::vector<bool>& values) {
        if (values.size() != size_) {
            resize(values.size());
        }
        data_ = values;
    }
    
    const std::vector<bool>& data() const { return data_; }
    std::vector<bool>& data() { return data_; }
    
    // BaseData 接口实现
    void resize(Index new_size) override {
        data_.resize(new_size, default_value_);
        size_ = new_size;
    }
    
    void clear() override {
        data_.clear();
        size_ = 0;
    }
    
    void reset() override {
        std::fill(data_.begin(), data_.end(), default_value_);
    }
    
    std::string type_name() const override {
        return "Bool";
    }
    
    std::unique_ptr<BaseData> clone() const override {
        auto cloned = std::make_unique<FieldData<bool>>(name_, location_, size_);
        cloned->data_ = data_;
        cloned->default_value_ = default_value_;
        return cloned;
    }
    
    bool default_value() const { return default_value_; }
    void set_default_value(bool value) { default_value_ = value; }

private:
    std::vector<bool> data_;
    bool default_value_;
    
    void check_index(Index index) const {
        if (index >= size_) {
            std::ostringstream oss;
            oss << "FieldData<bool> '" << name_ << "' index out of range: "
                << index << " >= " << size_;
            throw std::out_of_range(oss.str());
        }
    }
};

// ═══════════════════════════════════════════════════════════
// 类型别名（方便使用）
// ═══════════════════════════════════════════════════════════

using RealData = FieldData<Real>;              ///< 标量场
using IntData = FieldData<int>;                ///< 整型场
using BoolData = FieldData<bool>;              ///< 布尔场（特化版本）
using StringData = FieldData<std::string>;     ///< 字符串场
using VectorData = FieldData<Vector>;          ///< 矢量场
using TensorData = FieldData<DenseMatrix>;     ///< 张量场

} // namespace data
} // namespace fem
