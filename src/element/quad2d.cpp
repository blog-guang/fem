#include "element/quad2d.h"
#include <cmath>

namespace fem {

// 2x2 Gauss 积分点
static const Real GAUSS_1 = 1.0 / std::sqrt(3.0);  // ≈ 0.57735

static const QuadPoint QUAD_QUAD_POINTS[] = {
    {{-GAUSS_1, -GAUSS_1, 0.0}, 1.0},
    {{ GAUSS_1, -GAUSS_1, 0.0}, 1.0},
    {{ GAUSS_1,  GAUSS_1, 0.0}, 1.0},
    {{-GAUSS_1,  GAUSS_1, 0.0}, 1.0}
};

Span<const QuadPoint> Quad2D::quad_points() const {
    return {QUAD_QUAD_POINTS, 4};
}

void Quad2D::shape_functions(const Vec3& xi, Real* out) const {
    Real x = xi[0], y = xi[1];
    out[0] = 0.25 * (1.0 - x) * (1.0 - y);
    out[1] = 0.25 * (1.0 + x) * (1.0 - y);
    out[2] = 0.25 * (1.0 + x) * (1.0 + y);
    out[3] = 0.25 * (1.0 - x) * (1.0 + y);
}

void Quad2D::shape_gradients(const Vec3& xi, Real* out) const {
    Real x = xi[0], y = xi[1];
    //         dN/dxi              dN/deta          dN/dzeta
    out[0] = -0.25*(1.0-y);  out[1] = -0.25*(1.0-x);  out[2]  = 0.0;  // N0
    out[3] =  0.25*(1.0-y);  out[4] = -0.25*(1.0+x);  out[5]  = 0.0;  // N1
    out[6] =  0.25*(1.0+y);  out[7] =  0.25*(1.0+x);  out[8]  = 0.0;  // N2
    out[9] = -0.25*(1.0+y);  out[10]=  0.25*(1.0-x);  out[11] = 0.0;  // N3
}

}  // namespace fem
