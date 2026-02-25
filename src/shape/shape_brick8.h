#pragma once

#include "shape_function.h"

namespace fem {
namespace shape {

/**
 * Brick8ShapeFunction - 六面体8节点三线性单元
 * 
 * 自然坐标：(ξ, η, ζ) ∈ [-1,1]³
 * 
 * 形函数（三线性）：
 *   N_i = 1/8 * (1 + ξ*ξ_i) * (1 + η*η_i) * (1 + ζ*ζ_i)
 */
class Brick8ShapeFunction : public ShapeFunction3D {
public:
    void evaluate(const Vec3& xi, Vector& N) const override;
    void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const override;
    void getGaussPoints(int order, std::vector<Vec3>& points,
                       std::vector<Real>& weights) const override;
    
    int numNodes() const override { return 8; }
    ElementType elementType() const override { return ElementType::Brick8; }
};

}  // namespace shape
}  // namespace fem
