#include "physics/elasticity_material.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void ElasticityWithMaterial::initialize(std::size_t num_elements) {
    element_states_.clear();
    element_stresses_.clear();
    
    element_states_.reserve(num_elements);
    element_stresses_.reserve(num_elements);
    
    for (std::size_t i = 0; i < num_elements; ++i) {
        element_states_.push_back((*material_).createState());
        element_stresses_.push_back(Vector(3, 0.0));  // 2D: 3分量
    }
}

void ElasticityWithMaterial::compute_stiffness(Index elem_id, const Mesh& mesh,
                                              DenseMatrix& Ke) const {
    const Element& elem = mesh.element(elem_id);
    
    // 目前仅支持 Tri3
    if (elem.type() != ElementType::Tri3) {
        FEM_WARN("ElasticityWithMaterial: unsupported element type");
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    // 计算 B 矩阵
    DenseMatrix B(3, 6);
    Real area = compute_B_matrix(coords, B);
    
    if (area < 1e-15) {
        FEM_WARN("ElasticityWithMaterial: degenerate element");
        return;
    }
    
    // 获取材料切线刚度矩阵
    Vector strain_dummy(3, 0.0);  // 切线刚度可能依赖当前应变
    DenseMatrix D(3, 3);
    
    if (elem_id < element_states_.size()) {
        (*material_).computeTangent(strain_dummy, D, element_states_[elem_id]);
    } else {
        StateVariables temp_state = (*material_).createState();
        (*material_).computeTangent(strain_dummy, D, temp_state);
    }
    
    // Ke = B^T D B * area * thickness
    DenseMatrix DB = D * B;  // 3x6
    DenseMatrix Bt = B.transpose();  // 6x3
    DenseMatrix temp = Bt * DB;  // 6x6
    
    Ke = temp * (area * thickness_);
}

void ElasticityWithMaterial::update_stress(Index elem_id, const Mesh& mesh,
                                          const Vector& u_inc, Vector& Fe) {
    const Element& elem = mesh.element(elem_id);
    
    if (elem.type() != ElementType::Tri3) {
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    // 计算 B 矩阵
    DenseMatrix B(3, 6);
    Real area = compute_B_matrix(coords, B);
    
    // 提取单元位移增量
    Vector u_elem(6, 0.0);
    for (int i = 0; i < 3; ++i) {
        Index node_id = nodes[i];
        u_elem[2*i]     = u_inc[2*node_id];      // u_x
        u_elem[2*i + 1] = u_inc[2*node_id + 1];  // u_y
    }
    
    // 计算应变增量: Δε = B * Δu
    Vector strain_inc(3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            strain_inc[i] += B(i, j) * u_elem[j];
        }
    }
    
    // 更新应力
    if (elem_id >= element_stresses_.size()) {
        element_stresses_.resize(elem_id + 1, Vector(3, 0.0));
        element_states_.resize(elem_id + 1, (*material_).createState());
    }
    
    (*material_).computeStress(strain_inc, 
                              element_stresses_[elem_id], 
                              element_states_[elem_id]);
    
    // 计算内力: Fe = B^T * σ * area * thickness
    const Vector& stress = element_stresses_[elem_id];
    DenseMatrix Bt = B.transpose();  // 6x3
    
    Fe.resize(6, 0.0);
    for (std::size_t i = 0; i < 6; ++i) {
        Fe[i] = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
            Fe[i] += Bt(i, j) * stress[j];
        }
        Fe[i] *= area * thickness_;
    }
}

void ElasticityWithMaterial::compute_element(Index elem_id, const Mesh& mesh,
                                            const Vector& u_total,
                                            DenseMatrix& Ke, Vector& Fe) {
    // 1. 计算切线刚度
    compute_stiffness(elem_id, mesh, Ke);
    
    // 2. 计算内力（基于当前应力）
    const Element& elem = mesh.element(elem_id);
    const auto& nodes = elem.nodes();
    
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    DenseMatrix B(3, 6);
    Real area = compute_B_matrix(coords, B);
    
    if (elem_id < element_stresses_.size()) {
        const Vector& stress = element_stresses_[elem_id];
        DenseMatrix Bt = B.transpose();
        
        Fe.resize(6, 0.0);
        for (std::size_t i = 0; i < 6; ++i) {
            Fe[i] = 0.0;
            for (std::size_t j = 0; j < 3; ++j) {
                Fe[i] += Bt(i, j) * stress[j];
            }
            Fe[i] *= area * thickness_;
        }
    } else {
        Fe.resize(6, 0.0);
    }
}

Vector ElasticityWithMaterial::compute_stress(Index elem_id, const Mesh& mesh,
                                             const Vector& u_total) const {
    if (elem_id < element_stresses_.size()) {
        return element_stresses_[elem_id];
    }
    return Vector(3, 0.0);
}

const StateVariables& ElasticityWithMaterial::get_state(Index elem_id) const {
    static StateVariables dummy;
    if (elem_id < element_states_.size()) {
        return element_states_[elem_id];
    }
    return dummy;
}

void ElasticityWithMaterial::reset() {
    for (auto& state : element_states_) {
        state.reset();
    }
    for (auto& stress : element_stresses_) {
        stress.zero();
    }
}

Real ElasticityWithMaterial::compute_B_matrix(const Vec3 coords[3], 
                                             DenseMatrix& B) const {
    Real x1 = coords[0][0], y1 = coords[0][1];
    Real x2 = coords[1][0], y2 = coords[1][1];
    Real x3 = coords[2][0], y3 = coords[2][1];
    
    Real area = 0.5 * std::abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));
    
    if (area < 1e-15) {
        B.zero();
        return 0.0;
    }
    
    // B 矩阵（常应变三角形）
    Real inv_2A = 1.0 / (2.0 * area);
    
    Real b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    Real c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
    
    B.zero();
    
    // ε_xx = ∂u/∂x
    B(0, 0) = b1 * inv_2A;
    B(0, 2) = b2 * inv_2A;
    B(0, 4) = b3 * inv_2A;
    
    // ε_yy = ∂v/∂y
    B(1, 1) = c1 * inv_2A;
    B(1, 3) = c2 * inv_2A;
    B(1, 5) = c3 * inv_2A;
    
    // γ_xy = ∂u/∂y + ∂v/∂x
    B(2, 0) = c1 * inv_2A;
    B(2, 1) = b1 * inv_2A;
    B(2, 2) = c2 * inv_2A;
    B(2, 3) = b2 * inv_2A;
    B(2, 4) = c3 * inv_2A;
    B(2, 5) = b3 * inv_2A;
    
    return area;
}

}  // namespace physics
}  // namespace fem
