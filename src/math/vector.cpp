#include "math/vector.h"
#include "core/logger.h"
#include <algorithm>
#include <cmath>

namespace fem {

Vector Vector::operator-() const {
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = -data_[i];
    }
    return result;
}

Vector Vector::operator+(const Vector& other) const {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in addition");
    }
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_[i] + other[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in subtraction");
    }
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_[i] - other[i];
    }
    return result;
}

Vector Vector::operator*(Real scalar) const {
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_[i] * scalar;
    }
    return result;
}

Vector Vector::operator/(Real scalar) const {
    if (std::abs(scalar) < 1e-15) {
        throw std::runtime_error("Division by zero in vector");
    }
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_[i] / scalar;
    }
    return result;
}

Vector& Vector::operator+=(const Vector& other) {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in addition");
    }
    for (std::size_t i = 0; i < size(); ++i) {
        data_[i] += other[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& other) {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in subtraction");
    }
    for (std::size_t i = 0; i < size(); ++i) {
        data_[i] -= other[i];
    }
    return *this;
}

Vector& Vector::operator*=(Real scalar) {
    for (std::size_t i = 0; i < size(); ++i) {
        data_[i] *= scalar;
    }
    return *this;
}

Vector& Vector::operator/=(Real scalar) {
    if (std::abs(scalar) < 1e-15) {
        throw std::runtime_error("Division by zero in vector");
    }
    for (std::size_t i = 0; i < size(); ++i) {
        data_[i] /= scalar;
    }
    return *this;
}

Real Vector::dot(const Vector& other) const {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in dot product");
    }
    Real result = 0.0;
    for (std::size_t i = 0; i < size(); ++i) {
        result += data_[i] * other[i];
    }
    return result;
}

Vector Vector::hadamard(const Vector& other) const {
    if (size() != other.size()) {
        throw std::runtime_error("Vector size mismatch in Hadamard product");
    }
    Vector result(size());
    for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_[i] * other[i];
    }
    return result;
}

Real Vector::norm() const {
    return std::sqrt(norm_squared());
}

Real Vector::norm1() const {
    Real sum = 0.0;
    for (Real x : data_) {
        sum += std::abs(x);
    }
    return sum;
}

Real Vector::norm_inf() const {
    Real max_val = 0.0;
    for (Real x : data_) {
        max_val = std::max(max_val, std::abs(x));
    }
    return max_val;
}

Real Vector::norm_squared() const {
    Real sum = 0.0;
    for (Real x : data_) {
        sum += x * x;
    }
    return sum;
}

void Vector::normalize() {
    Real n = norm();
    if (n < 1e-15) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    *this /= n;
}

Vector Vector::normalized() const {
    Vector result = *this;
    result.normalize();
    return result;
}

void Vector::print(const std::string& name) const {
    if (!name.empty()) {
        FEM_INFO("Vector " + name + " (size=" + std::to_string(size()) + "):");
    } else {
        FEM_INFO("Vector (size=" + std::to_string(size()) + "):");
    }
    
    std::size_t print_limit = 10;
    for (std::size_t i = 0; i < std::min(size(), print_limit); ++i) {
        FEM_INFO("  [" + std::to_string(i) + "] = " + fmt_sci(data_[i]));
    }
    if (size() > print_limit) {
        FEM_INFO("  ... (" + std::to_string(size() - print_limit) + " more)");
    }
}

}  // namespace fem
