#include "shape_function_factory.h"
#include <stdexcept>

namespace fem {
namespace shape {

ShapeFunctionPtr ShapeFunctionFactory::create(ElementType type) {
    switch (type) {
        case ElementType::Tri3:
            return std::make_shared<Tri3ShapeFunction>();
        
        case ElementType::Quad4:
            return std::make_shared<Quad4ShapeFunction>();
        
        case ElementType::Tet4:
            return std::make_shared<Tet4ShapeFunction>();
        
        case ElementType::Brick8:
            return std::make_shared<Brick8ShapeFunction>();
        
        // 未实现的单元类型
        default:
            throw std::invalid_argument("Unsupported element type for shape function");
    }
}

void ShapeFunctionFactory::evaluate(ElementType type, const Vec3& xi, Vector& N) {
    auto shape_func = create(type);
    shape_func->evaluate(xi, N);
}

void ShapeFunctionFactory::evaluateDerivatives(ElementType type, const Vec3& xi, 
                                              DenseMatrix& dN) {
    auto shape_func = create(type);
    shape_func->evaluateDerivatives(xi, dN);
}

void ShapeFunctionFactory::getGaussPoints(ElementType type, int order,
                                         std::vector<Vec3>& points,
                                         std::vector<Real>& weights) {
    auto shape_func = create(type);
    shape_func->getGaussPoints(order, points, weights);
}

}  // namespace shape
}  // namespace fem
