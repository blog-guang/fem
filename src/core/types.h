#pragma once

#include <cstddef>
#include <cstdint>
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

namespace fem {

// ── 基础类型 ──
using Real   = double;
using Index  = std::size_t;
using Vec2   = std::array<Real, 2>;
using Vec3   = std::array<Real, 3>;
using Vector = std::vector<Real>;

// ── 单元节点上限 (Hex=8, 覆盖所有常用单元) ──
constexpr std::size_t MAX_NODES = 8;

// ── 枚举 ──
enum class ElementType : uint8_t {
    // Legacy (Phase 1/2 compatibility)
    Triangle2D,      // 三角形 (3 节点)
    Quad2D,          // 四边形 (4 节点)
    Tetrahedron3D,   // 四面体 (4 节点)
    Hexahedron3D,    // 六面体 (8 节点)
    
    // Mesh V2 (新架构)
    Node,            // 0D
    Edge2, Edge3,    // 1D
    Tri3, Tri6,      // 2D 三角形
    Quad4, Quad8,    // 2D 四边形
    Tet4, Tet10,     // 3D 四面体
    Brick8, Brick20  // 3D 六面体
};

enum class FieldType : uint8_t {
    Scalar,          // 标量场
    Vector2D,        // 2D 矢量场
    Vector3D         // 3D 矢量场
};

enum class BCType : uint8_t {
    Dirichlet,       // 第一类: 固定值
    Neumann          // 第二类: 固定梯度/通量
};

// ── 简单 Span (替代 C++20 std::span, GCC<20 兼容) ──
template<typename T>
struct Span {
    const T*    data{nullptr};
    std::size_t size{0};

    const T* begin() const { return data; }
    const T* end()   const { return data + size; }
    const T& operator[](std::size_t i) const { return data[i]; }
};

// ── 格式化辅助 (GCC<13 兼容, 替代 std::format) ──
inline std::string fmt_sci(Real val, int prec = 2) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(prec) << val;
    return oss.str();
}

}  // namespace fem
