#include "physics/heat_unified.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

void HeatConductionUnified::compute_element(Index elem_id, const Mesh& mesh,
                                           DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    ElementType elem_type = elem.type();
    
    // 获取节点坐标
    std::vector<Vec3> coords = get_element_coords(elem_id, mesh);
    Index n_nodes = coords.size();
    
    // 创建形函数对象
    auto shape_func = shape::ShapeFunctionFactory::create(elem_type);
    if (!shape_func) {
        FEM_WARN("HeatConductionUnified: unsupported element type");
        return;
    }
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = get_gauss_order(elem_type);
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    // 初始化单元矩阵
    Ke.resize(n_nodes, n_nodes);
    Ke.zero();
    Fe.resize(n_nodes);
    for (Index i = 0; i < n_nodes; ++i) {
        Fe[i] = 0.0;
    }
    
    // 高斯积分循环
    for (size_t gp = 0; gp < gauss_points.size(); ++gp) {
        const Vec3& xi = gauss_points[gp];
        Real w = weights[gp];
        
        // 计算形函数值和导数
        Vector N;
        shape_func->evaluate(xi, N);
        
        DenseMatrix dN_dxyz;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dxyz);
        
        // 计算雅可比行列式
        DenseMatrix J = shape_func->computeJacobian(xi, coords);
        Real det_J = compute_jacobian_determinant(J);
        
        if (std::abs(det_J) < 1e-15) {
            FEM_WARN("HeatConductionUnified: degenerate element (det(J) near zero)");
            continue;
        }
        
        Real dV = w * std::abs(det_J);  // dΩ = |det(J)| * w
        
        // 装配刚度矩阵: Ke += k * ∇N^T * ∇N * dV
        int dim = shape_func->dimension();
        for (Index i = 0; i < n_nodes; ++i) {
            for (Index j = 0; j < n_nodes; ++j) {
                Real grad_dot = 0.0;
                for (int d = 0; d < dim; ++d) {
                    grad_dot += dN_dxyz(i, d) * dN_dxyz(j, d);
                }
                Ke(i, j) += k_ * grad_dot * dV;
            }
        }
        
        // 装配载荷向量: Fe += Q * N * dV
        for (Index i = 0; i < n_nodes; ++i) {
            Fe[i] += Q_ * N[i] * dV;
        }
    }
}

void HeatConductionUnified::compute_stiffness(Index elem_id, const Mesh& mesh,
                                             DenseMatrix& Ke) const {
    Vector Fe_dummy;
    compute_element(elem_id, mesh, Ke, Fe_dummy);
}

void HeatConductionUnified::compute_load(Index elem_id, const Mesh& mesh,
                                        Vector& Fe) const {
    DenseMatrix Ke_dummy;
    compute_element(elem_id, mesh, Ke_dummy, Fe);
}

} // namespace physics
} // namespace fem
