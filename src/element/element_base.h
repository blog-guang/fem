#pragma once

#include "core/types.h"
#include <memory>

namespace fem {

// ── 积分点 (POD) ──
struct QuadPoint {
    Vec3 xi;        // 参考坐标
    Real weight;    // 权重
};

// ── 有限元基类 ──
// 接口薄: 只暴露形函数 + 积分点。具体单元实现在子类。
class ElementBase {
public:
    virtual ~ElementBase() = default;

    virtual uint8_t num_nodes() const = 0;

    // 积分点集 (静态数据, 无分配)
    virtual Span<const QuadPoint> quad_points() const = 0;

    // 形函数值 N_i(xi), 输出写入调用方提供的缓冲区
    virtual void shape_functions(const Vec3& xi, Real* out) const = 0;

    // 形函数梯度 dN_i/dxi, 输出写入调用方提供的缓冲区
    // out[i*3 + d] = dN_i / d(xi_d)
    virtual void shape_gradients(const Vec3& xi, Real* out) const = 0;
};

// ── 工厂 ──
[[nodiscard]] std::unique_ptr<ElementBase> create_element(ElementType type);

}  // namespace fem
