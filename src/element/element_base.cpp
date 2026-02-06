#include "element/element_base.h"
#include "element/triangle2d.h"
#include "element/quad2d.h"
#include <stdexcept>

namespace fem {

std::unique_ptr<ElementBase> create_element(ElementType type) {
    switch (type) {
        case ElementType::Triangle2D:   return std::make_unique<Triangle2D>();
        case ElementType::Quad2D:       return std::make_unique<Quad2D>();
        default:
            throw std::invalid_argument(
                "create_element: type " + std::to_string(static_cast<int>(type)) + " not implemented");
    }
}

}  // namespace fem
