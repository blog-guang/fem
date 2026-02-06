#pragma once

#include "core/types.h"
#include <vector>
#include <cmath>
#include <stdexcept>

namespace fem {

/**
 * Vector: 动态大小的向量类
 * 
 * 支持基本的向量运算：加减、数乘、点积、范数等
 */
class Vector {
public:
    // ═══ 构造 ═══
    Vector() = default;
    explicit Vector(std::size_t size, Real value = 0.0) 
        : data_(size, value) {}
    
    Vector(const std::vector<Real>& data) : data_(data) {}
    Vector(std::vector<Real>&& data) : data_(std::move(data)) {}
    
    // ═══ 大小与访问 ═══
    std::size_t size() const { return data_.size(); }
    void resize(std::size_t n, Real value = 0.0) { data_.resize(n, value); }
    void clear() { data_.clear(); }
    
    Real& operator[](std::size_t i) { return data_[i]; }
    const Real& operator[](std::size_t i) const { return data_[i]; }
    
    Real& at(std::size_t i) { 
        if (i >= data_.size()) throw std::out_of_range("Vector index out of range");
        return data_[i]; 
    }
    const Real& at(std::size_t i) const { 
        if (i >= data_.size()) throw std::out_of_range("Vector index out of range");
        return data_[i]; 
    }
    
    Real* data() { return data_.data(); }
    const Real* data() const { return data_.data(); }
    
    std::vector<Real>& raw() { return data_; }
    const std::vector<Real>& raw() const { return data_; }
    
    // ═══ 赋值 ═══
    void fill(Real value) { std::fill(data_.begin(), data_.end(), value); }
    void zero() { fill(0.0); }
    
    // ═══ 向量运算 ═══
    // v + w
    Vector operator+(const Vector& other) const;
    // v - w
    Vector operator-(const Vector& other) const;
    // a * v
    Vector operator*(Real scalar) const;
    // v / a
    Vector operator/(Real scalar) const;
    
    // v += w
    Vector& operator+=(const Vector& other);
    // v -= w
    Vector& operator-=(const Vector& other);
    // v *= a
    Vector& operator*=(Real scalar);
    // v /= a
    Vector& operator/=(Real scalar);
    
    // 点积: v · w
    Real dot(const Vector& other) const;
    
    // 范数
    Real norm() const;           // L2 范数
    Real norm1() const;          // L1 范数
    Real norm_inf() const;       // L∞ 范数
    Real norm_squared() const;   // ||v||^2
    
    // 归一化
    void normalize();
    Vector normalized() const;
    
    // ═══ 工具函数 ═══
    void print(const std::string& name = "") const;
    
private:
    std::vector<Real> data_;
};

// ═══ 全局运算符 ═══
// a * v (标量在左)
inline Vector operator*(Real scalar, const Vector& v) {
    return v * scalar;
}

}  // namespace fem
