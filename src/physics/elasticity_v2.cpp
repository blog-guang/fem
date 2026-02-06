#include "physics/elasticity_v2.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void Elasticity2D::compute_D_matrix() {
    D_ = DenseMatrix(3, 3, 0.0);
    
    if (plane_type_ == PlaneType::PlaneStress) {
        // 平面应力
        Real factor = E_ / (1.0 - nu_ * nu_);
        D_(0, 0) = factor;
        D_(0, 1) = factor * nu_;
        D_(1, 0) = factor * nu_;
        D_(1, 1) = factor;
        D_(2, 2) = factor * (1.0 - nu_) / 2.0;
    } else {
        // 平面应变
        Real factor = E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
        D_(0, 0) = factor * (1.0 - nu_);
        D_(0, 1) = factor * nu_;
        D_(1, 0) = factor * nu_;
        D_(1, 1) = factor * (1.0 - nu_);
        D_(2, 2) = factor * (1.0 - 2.0 * nu_) / 2.0;
    }
}

void Elasticity2D::compute_element(Index elem_id, const Mesh& mesh,
                                  DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    
    // 目前仅支持 Tri3
    if (elem.type() != ElementType::Tri3) {
        FEM_WARN("Elasticity2D: unsupported element type, skipping");
        return;
    }
    
    const auto& nodes = elem.nodes();
    Vec3 coords[3];
    for (int i = 0; i < 3; ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    // 计算 B 矩阵和面积
    DenseMatrix B(3, 6);  // 3 (ε分量) x 6 (u_x1, u_y1, u_x2, u_y2, u_x3, u_y3)
    Real area = compute_tri3_B_matrix(coords, B);
    
    if (area < 1e-15) {
        FEM_WARN("Elasticity2D: degenerate element (area near zero)");
        return;
    }
    
    // 计算刚度矩阵: Ke = ∫_Ω B^T D B dΩ = B^T D B * area
    DenseMatrix DB = D_ * B;  // 3x6
    DenseMatrix Bt = B.transpose();  // 6x3
    
    Ke = Bt * DB * area;  // 6x6
    
    // 载荷向量 (暂不考虑体力)
    for (std::size_t i = 0; i < Fe.size(); ++i) {
        Fe[i] = 0.0;
    }
}

void Elasticity2D::compute_stiffness(Index elem_id, const Mesh& mesh,
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
    
    DenseMatrix B(3, 6);
    Real area = compute_tri3_B_matrix(coords, B);
    
    if (area < 1e-15) return;
    
    DenseMatrix DB = D_ * B;
    DenseMatrix Bt = B.transpose();
    
    Ke = Bt * DB * area;
}

Real Elasticity2D::compute_tri3_B_matrix(const Vec3* coords, DenseMatrix& B) const {
    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];
    
    // 计算面积
    Real detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    Real area = 0.5 * std::abs(detJ);
    
    if (area < 1e-15) {
        return area;
    }
    
    // 计算形函数梯度
    Real dN_dx[3], dN_dy[3];
    
    // ∇N0
    dN_dx[0] = (y1 - y2) / (2.0 * area);
    dN_dy[0] = (x2 - x1) / (2.0 * area);
    
    // ∇N1
    dN_dx[1] = (y2 - y0) / (2.0 * area);
    dN_dy[1] = (x0 - x2) / (2.0 * area);
    
    // ∇N2
    dN_dx[2] = (y0 - y1) / (2.0 * area);
    dN_dy[2] = (x1 - x0) / (2.0 * area);
    
    // 构造 B 矩阵 (3x6)
    // B = [[dN1/dx,  0,       dN2/dx, 0,       dN3/dx, 0      ],
    //      [0,       dN1/dy,  0,      dN2/dy,  0,      dN3/dy ],
    //      [dN1/dy,  dN1/dx,  dN2/dy, dN2/dx,  dN3/dy, dN3/dx]]
    
    B.zero();
    
    for (int i = 0; i < 3; ++i) {
        int col_u = i * 2;      // u_x 列
        int col_v = i * 2 + 1;  // u_y 列
        
        B(0, col_u) = dN_dx[i];  // ε_xx = ∂u_x/∂x
        B(1, col_v) = dN_dy[i];  // ε_yy = ∂u_y/∂y
        B(2, col_u) = dN_dy[i];  // γ_xy = ∂u_x/∂y + ∂u_y/∂x
        B(2, col_v) = dN_dx[i];
    }
    
    return area;
}

} // namespace physics
} // namespace fem
