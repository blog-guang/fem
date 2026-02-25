#pragma once

#include "shape_function.h"

namespace fem {
namespace shape {

/**
 * Quad4ShapeFunction - 四边形4节点双线性单元
 * 
 * 节点编号：
 *   3---2
 *   |   |
 *   |   |
 *   0---1
 * 
 * 自然坐标：(ξ, η) ∈ [-1,1] × [-1,1]
 * 
 * 形函数（双线性）：
 *   N_i = 1/4 * (1 + ξ*ξ_i) * (1 + η*η_i)
 *   其中 (ξ_i, η_i) 是节点i的自然坐标
 */
class Quad4ShapeFunction : public ShapeFunction2D {
public:
    void evaluate(const Vec3& xi, Vector& N) const override;
    
    void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const override;
    
    void getGaussPoints(int order,
                       std::vector<Vec3>& points,
                       std::vector<Real>& weights) const override;
    
    int numNodes() const override { return 4; }
    ElementType elementType() const override { return ElementType::Quad4; }
};

}  // namespace shape
}  // namespace fem
