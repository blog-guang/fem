#include "physics/heat.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void HeatConduction::compute_element(Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    
    // 目前仅支持 Tri3
    if (elem.type() != ElementType::Tri3) {
        FEM_WARN("HeatConduction: unsupported element type, skipping");
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    // 计算形函数梯度和面积
    Real grad[3][2];
    Real area = compute_tri3_gradients(coords, grad);
    
    if (area < 1e-15) {
        FEM_WARN("HeatConduction: degenerate element (area near zero)");
        return;
    }
    
    // 计算刚度矩阵: Ke[i][j] = k * ∫_Ω ∇Ni · ∇Nj dΩ
    //                         = k * (∇Ni · ∇Nj) * area
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke(i, j) = k_ * (grad[i][0]*grad[j][0] + grad[i][1]*grad[j][1]) * area;
        }
    }
    
    // 计算载荷向量: Fe[i] = Q * ∫_Ω Ni dΩ
    //                      = Q * area / 3.0  (Ni 在单元内平均值为 1/3)
    for (int i = 0; i < 3; ++i) {
        Fe[i] = Q_ * area / 3.0;
    }
}

void HeatConduction::compute_stiffness(Index elem_id, const Mesh& mesh,
                                      DenseMatrix& Ke) const {
    const Element& elem = mesh.element(elem_id);
    
    if (elem.type() != ElementType::Tri3) {
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    Real grad[3][2];
    Real area = compute_tri3_gradients(coords, grad);
    
    if (area < 1e-15) return;
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke(i, j) = k_ * (grad[i][0]*grad[j][0] + grad[i][1]*grad[j][1]) * area;
        }
    }
}

void HeatConduction::compute_load(Index elem_id, const Mesh& mesh,
                                 Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    
    if (elem.type() != ElementType::Tri3) {
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    Real grad[3][2];  // 占位
    Real area = compute_tri3_gradients(coords, grad);
    
    if (area < 1e-15) return;
    
    for (int i = 0; i < 3; ++i) {
        Fe[i] = Q_ * area / 3.0;
    }
}

Real HeatConduction::compute_tri3_gradients(const Vec3* coords, Real grad[][2]) const {
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];
    
    // 计算雅可比行列式和面积
    Real detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    Real area = 0.5 * std::abs(detJ);
    
    if (area < 1e-15) {
        return area;
    }
    
    // 计算形函数梯度
    // ∇N0
    grad[0][0] = (y1 - y2) / (2.0 * area);
    grad[0][1] = (x2 - x1) / (2.0 * area);
    
    // ∇N1
    grad[1][0] = (y2 - y0) / (2.0 * area);
    grad[1][1] = (x0 - x2) / (2.0 * area);
    
    // ∇N2
    grad[2][0] = (y0 - y1) / (2.0 * area);
    grad[2][1] = (x1 - x0) / (2.0 * area);
    
    return area;
}

} // namespace physics
} // namespace fem
