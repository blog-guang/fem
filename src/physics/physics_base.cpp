#include "physics/physics_base.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

std::vector<Vec3> PhysicsBase::get_element_coords(Index elem_id, const Mesh& mesh) const {
    const Element& elem = mesh.element(elem_id);
    const auto& nodes = elem.nodes();
    
    std::vector<Vec3> coords(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        coords[i] = mesh.node(nodes[i]).coords();
    }
    
    return coords;
}

bool PhysicsBase::check_dimension(const Element& elem, int expected_dim) const {
    ElementType type = elem.type();
    
    // 2D 单元
    if (expected_dim == 2) {
        return type == ElementType::Tri3 || type == ElementType::Quad4;
    }
    
    // 3D 单元
    if (expected_dim == 3) {
        return type == ElementType::Tet4 || type == ElementType::Brick8;
    }
    
    return false;
}

std::tuple<Index, Index, int> PhysicsBase::get_element_info(
    const Element& elem, int dofs_per_node) const {
    
    Index n_nodes = elem.nodes().size();
    Index n_dofs = n_nodes * dofs_per_node;
    
    // 推断维度
    int dim = 0;
    ElementType type = elem.type();
    if (type == ElementType::Tri3 || type == ElementType::Quad4) {
        dim = 2;
    } else if (type == ElementType::Tet4 || type == ElementType::Brick8) {
        dim = 3;
    }
    
    return {n_nodes, n_dofs, dim};
}

int PhysicsBase::get_gauss_order(ElementType elem_type) const {
    // 线性单元用2阶积分（精确）
    // 未来可扩展为高阶单元用更高阶积分
    return 2;
}

DenseMatrix PhysicsBase::build_B_matrix_2D(const DenseMatrix& dN_dxyz) const {
    Index n_nodes = dN_dxyz.rows();
    Index n_dofs = n_nodes * 2;
    
    DenseMatrix B(3, n_dofs, 0.0);
    
    for (Index i = 0; i < n_nodes; ++i) {
        Real dNdx = dN_dxyz(i, 0);
        Real dNdy = dN_dxyz(i, 1);
        
        Index col_x = i * 2;
        Index col_y = i * 2 + 1;
        
        B(0, col_x) = dNdx;      // ε_xx = ∂u_x/∂x
        B(1, col_y) = dNdy;      // ε_yy = ∂u_y/∂y
        B(2, col_x) = dNdy;      // γ_xy = ∂u_x/∂y + ∂u_y/∂x
        B(2, col_y) = dNdx;
    }
    
    return B;
}

DenseMatrix PhysicsBase::build_B_matrix_3D(const DenseMatrix& dN_dxyz) const {
    Index n_nodes = dN_dxyz.rows();
    Index n_dofs = n_nodes * 3;
    
    DenseMatrix B(6, n_dofs, 0.0);
    
    for (Index i = 0; i < n_nodes; ++i) {
        Real dNdx = dN_dxyz(i, 0);
        Real dNdy = dN_dxyz(i, 1);
        Real dNdz = dN_dxyz(i, 2);
        
        Index col_x = i * 3;
        Index col_y = i * 3 + 1;
        Index col_z = i * 3 + 2;
        
        B(0, col_x) = dNdx;      // ε_xx = ∂u_x/∂x
        B(1, col_y) = dNdy;      // ε_yy = ∂u_y/∂y
        B(2, col_z) = dNdz;      // ε_zz = ∂u_z/∂z
        B(3, col_y) = dNdz;      // γ_yz = ∂u_y/∂z + ∂u_z/∂y
        B(3, col_z) = dNdy;
        B(4, col_x) = dNdz;      // γ_xz = ∂u_x/∂z + ∂u_z/∂x
        B(4, col_z) = dNdx;
        B(5, col_x) = dNdy;      // γ_xy = ∂u_x/∂y + ∂u_y/∂x
        B(5, col_y) = dNdx;
    }
    
    return B;
}

Real PhysicsBase::compute_jacobian_determinant(const DenseMatrix& J) const {
    Index dim = J.rows();
    
    if (dim == 2) {
        // 2D: det(J) = J00*J11 - J01*J10
        return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    } else if (dim == 3) {
        // 3D: det(J) = 萨吕斯法则
        return J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1))
             - J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0))
             + J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));
    }
    
    return 0.0;
}

} // namespace physics
} // namespace fem
