#pragma once

#include "element/element_base.h"

namespace fem {

// ── 双线性四边形单元 (Q1, 2D, 4 节点) ──
//
// 参考单元 [-1,1]x[-1,1]:
//   (-1,1) ── (1,1)
//     |        |
//   (-1,-1)── (1,-1)
//
// 节点编号: 0=(-1,-1), 1=(1,-1), 2=(1,1), 3=(-1,1)  (逆时针)
//
// 形函数:
//   N0 = 0.25*(1-xi)*(1-eta)
//   N1 = 0.25*(1+xi)*(1-eta)
//   N2 = 0.25*(1+xi)*(1+eta)
//   N3 = 0.25*(1-xi)*(1+eta)
//
class Quad2D : public ElementBase {
public:
    uint8_t num_nodes() const override { return 4; }

    Span<const QuadPoint> quad_points() const override;
    void shape_functions(const Vec3& xi, Real* out) const override;
    void shape_gradients(const Vec3& xi, Real* out) const override;
};

}  // namespace fem
