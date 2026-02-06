#include "element/triangle2d.h"

namespace fem {

// 1 点 Gauss 积分, 权重 = 参考三角形面积 0.5
static const QuadPoint TRI_QUAD_POINTS[] = {
    {{1.0/3.0, 1.0/3.0, 0.0}, 0.5}
};

Span<const QuadPoint> Triangle2D::quad_points() const {
    return {TRI_QUAD_POINTS, 1};
}

void Triangle2D::shape_functions(const Vec3& xi, Real* out) const {
    // N0 = 1 - xi - eta
    // N1 = xi
    // N2 = eta
    out[0] = 1.0 - xi[0] - xi[1];
    out[1] = xi[0];
    out[2] = xi[1];
}

void Triangle2D::shape_gradients(const Vec3& /*xi*/, Real* out) const {
    // 梯度为常数:
    //         dN/dxi   dN/deta   dN/dzeta
    // N0:     -1       -1        0
    // N1:      1        0        0
    // N2:      0        1        0
    // out[i*3 + d]
    out[0] = -1.0; out[1] = -1.0; out[2] = 0.0;   // N0
    out[3] =  1.0; out[4] =  0.0; out[5] = 0.0;   // N1
    out[6] =  0.0; out[7] =  1.0; out[8] = 0.0;   // N2
}

}  // namespace fem
