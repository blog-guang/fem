#pragma once

#include "shape_function.h"

namespace fem {
namespace shape {

/**
 * Tet4ShapeFunction - 四面体4节点线性单元
 * 
 * 自然坐标：(ξ, η, ζ) ∈ [0,1], ξ+η+ζ ≤ 1
 * 
 * 形函数：
 *   N1 = 1 - ξ - η - ζ
 *   N2 = ξ
 *   N3 = η
 *   N4 = ζ
 */
class Tet4ShapeFunction : public ShapeFunction3D {
public:
    void evaluate(const Vec3& xi, Vector& N) const override;
    void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const override;
    void getGaussPoints(int order, std::vector<Vec3>& points,
                       std::vector<Real>& weights) const override;
    
    int numNodes() const override { return 4; }
    ElementType elementType() const override { return ElementType::Tet4; }
};

}  // namespace shape
}  // namespace fem
