#include "state_variables.h"
#include <iostream>
#include <iomanip>

namespace fem {
namespace constitutive {

StateVariables::StateVariables(std::size_t strain_size)
    : plastic_strain(strain_size, 0.0),
      back_stress(strain_size, 0.0),
      equiv_plastic_strain(0.0),
      damage(0.0)
{}

void StateVariables::setScalar(const std::string& name, Real value) {
    scalar_vars[name] = value;
}

Real StateVariables::getScalar(const std::string& name, Real default_val) const {
    auto it = scalar_vars.find(name);
    return (it != scalar_vars.end()) ? it->second : default_val;
}

void StateVariables::setTensor(const std::string& name, const Vector& value) {
    tensor_vars[name] = value;
}

Vector StateVariables::getTensor(const std::string& name) const {
    auto it = tensor_vars.find(name);
    if (it != tensor_vars.end()) {
        return it->second;
    }
    return Vector();  // 空向量
}

bool StateVariables::hasTensor(const std::string& name) const {
    return tensor_vars.find(name) != tensor_vars.end();
}

void StateVariables::reset() {
    plastic_strain.zero();
    back_stress.zero();
    equiv_plastic_strain = 0.0;
    damage = 0.0;
    scalar_vars.clear();
    tensor_vars.clear();
}

StateVariables StateVariables::clone() const {
    StateVariables copy;
    copy.plastic_strain = plastic_strain;
    copy.equiv_plastic_strain = equiv_plastic_strain;
    copy.damage = damage;
    copy.back_stress = back_stress;
    copy.scalar_vars = scalar_vars;
    copy.tensor_vars = tensor_vars;
    return copy;
}

void StateVariables::serialize(std::ostream& os) const {
    // 简单二进制序列化 (生产环境可用HDF5/Protobuf)
    auto write_vector = [&](const Vector& v) {
        std::size_t n = v.size();
        os.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n > 0) {
            os.write(reinterpret_cast<const char*>(v.data()), n * sizeof(Real));
        }
    };
    
    // 预定义成员
    write_vector(plastic_strain);
    os.write(reinterpret_cast<const char*>(&equiv_plastic_strain), sizeof(Real));
    os.write(reinterpret_cast<const char*>(&damage), sizeof(Real));
    write_vector(back_stress);
    
    // 扩展标量
    std::size_t n_scalars = scalar_vars.size();
    os.write(reinterpret_cast<const char*>(&n_scalars), sizeof(n_scalars));
    for (const auto& [name, val] : scalar_vars) {
        std::size_t name_len = name.size();
        os.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
        os.write(name.data(), name_len);
        os.write(reinterpret_cast<const char*>(&val), sizeof(Real));
    }
    
    // 扩展张量
    std::size_t n_tensors = tensor_vars.size();
    os.write(reinterpret_cast<const char*>(&n_tensors), sizeof(n_tensors));
    for (const auto& [name, vec] : tensor_vars) {
        std::size_t name_len = name.size();
        os.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
        os.write(name.data(), name_len);
        write_vector(vec);
    }
}

void StateVariables::deserialize(std::istream& is) {
    auto read_vector = [&](Vector& v) {
        std::size_t n;
        is.read(reinterpret_cast<char*>(&n), sizeof(n));
        v.resize(n);
        if (n > 0) {
            is.read(reinterpret_cast<char*>(v.data()), n * sizeof(Real));
        }
    };
    
    // 预定义成员
    read_vector(plastic_strain);
    is.read(reinterpret_cast<char*>(&equiv_plastic_strain), sizeof(Real));
    is.read(reinterpret_cast<char*>(&damage), sizeof(Real));
    read_vector(back_stress);
    
    // 扩展标量
    std::size_t n_scalars;
    is.read(reinterpret_cast<char*>(&n_scalars), sizeof(n_scalars));
    scalar_vars.clear();
    for (std::size_t i = 0; i < n_scalars; ++i) {
        std::size_t name_len;
        is.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
        std::string name(name_len, '\0');
        is.read(&name[0], name_len);
        Real val;
        is.read(reinterpret_cast<char*>(&val), sizeof(Real));
        scalar_vars[name] = val;
    }
    
    // 扩展张量
    std::size_t n_tensors;
    is.read(reinterpret_cast<char*>(&n_tensors), sizeof(n_tensors));
    tensor_vars.clear();
    for (std::size_t i = 0; i < n_tensors; ++i) {
        std::size_t name_len;
        is.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
        std::string name(name_len, '\0');
        is.read(&name[0], name_len);
        Vector vec;
        read_vector(vec);
        tensor_vars[name] = vec;
    }
}

void StateVariables::print(const std::string& prefix) const {
    std::cout << prefix << "StateVariables:\n";
    std::cout << prefix << "  equiv_plastic_strain = " << equiv_plastic_strain << "\n";
    std::cout << prefix << "  damage = " << damage << "\n";
    
    if (plastic_strain.size() > 0) {
        std::cout << prefix << "  plastic_strain = [";
        for (std::size_t i = 0; i < plastic_strain.size(); ++i) {
            std::cout << plastic_strain[i];
            if (i < plastic_strain.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    if (!scalar_vars.empty()) {
        std::cout << prefix << "  scalar_vars:\n";
        for (const auto& [name, val] : scalar_vars) {
            std::cout << prefix << "    " << name << " = " << val << "\n";
        }
    }
}

}  // namespace constitutive
}  // namespace fem
