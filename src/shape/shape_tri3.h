#pragma once

#include "shape_function.h"

namespace fem {
namespace shape {

/**
 * Tri3ShapeFunction - 三角形3节点线性单元
 * 
 * 节点编号：
 *   2
 *   |\
 *   | \
 *   |  \
 *   0---1
 * 
 * 自然坐标：(ξ, η) ∈ [0,1], ξ+η ≤ 1
 * 
 * 形函数：
 *   N1 = 1 - ξ - η
 *   N2 = ξ
 *   N3 = η
 */
class Tri3ShapeFunction : public ShapeFunction2D {
public:
    void evaluate(const Vec3& xi, Vector& N) const override;
    
    void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const override;
    
    void getGaussPoints(int order,
                       std::vector<Vec3>& points,
                       std::vector<Real>& weights) const override;
    
    int numNodes() const override { return 3; }
    ElementType elementType() const override { return ElementType::Tri3; }
};

}  // namespace shape
}  // namespace fem
