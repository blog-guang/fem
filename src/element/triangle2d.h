#pragma once

#include "element/element_base.h"

namespace fem {

// ── 线性三角形单元 (P1, 2D, 3 节点) ──
//
// 参考单元:
//   (0,1)
//     |`.
//     |  `.
//     |    `.
//   (0,0)──(1,0)
//
// 形函数: N0 = 1-xi-eta, N1 = xi, N2 = eta
// 梯度为常数 (线性单元)
//
class Triangle2D : public ElementBase {
public:
    uint8_t num_nodes() const override { return 3; }

    Span<const QuadPoint> quad_points() const override;
    void shape_functions(const Vec3& xi, Real* out) const override;
    void shape_gradients(const Vec3& xi, Real* out) const override;
};

}  // namespace fem
