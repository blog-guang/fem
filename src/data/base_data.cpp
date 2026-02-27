#include "data/base_data.h"
#include <iostream>

namespace fem {
namespace data {

void BaseData::print_info() const {
    std::cout << "Field: " << name_ << "\n";
    std::cout << "  Type: " << type_name() << "\n";
    std::cout << "  Location: " << to_string(location_) << "\n";
    std::cout << "  Size: " << size_ << "\n";
}

} // namespace data
} // namespace fem
