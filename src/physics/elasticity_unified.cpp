#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace physics {

// ═══════════════════════════════════════════════════════════
// 构造函数实现
// ═══════════════════════════════════════════════════════════

ElasticityUnified::ElasticityUnified(Real youngs_modulus, Real poissons_ratio,
                                     PlaneType plane_type)
    : own_material_(true), dimension_(2), plane_type_(plane_type) {
    // 内部创建 2D 各向同性弹性材料
    bool plane_stress = (plane_type == PlaneType::PlaneStress);
    material_ = new constitutive::IsotropicElastic(
        youngs_modulus, 
        poissons_ratio, 
        2,  // dimension
        plane_stress
    );
}

ElasticityUnified::ElasticityUnified(Real youngs_modulus, Real poissons_ratio, bool use_3d)
    : own_material_(true), dimension_(use_3d ? 3 : 2), plane_type_(PlaneType::PlaneStress) {
    // 内部创建 2D 或 3D 各向同性弹性材料
    if (use_3d) {
        material_ = new constitutive::IsotropicElastic(
            youngs_modulus, 
            poissons_ratio, 
            3  // dimension
        );
    } else {
        material_ = new constitutive::IsotropicElastic(
            youngs_modulus, 
            poissons_ratio, 
            2,     // dimension
            true   // plane_stress
        );
    }
}

ElasticityUnified::ElasticityUnified(constitutive::Material* material, int dimension)
    : material_(material), own_material_(false), dimension_(dimension), 
      plane_type_(PlaneType::PlaneStress) {
    if (!material_) {
        FEM_ERROR("ElasticityUnified: material pointer is null");
    }
    if (dimension_ != 2 && dimension_ != 3) {
        FEM_ERROR("ElasticityUnified: dimension must be 2 or 3");
    }
}

ElasticityUnified::~ElasticityUnified() {
    if (own_material_ && material_) {
        delete material_;
        material_ = nullptr;
    }
}

// ═══════════════════════════════════════════════════════════
// 核心计算函数
// ═══════════════════════════════════════════════════════════

void ElasticityUnified::compute_element(Index elem_id, const Mesh& mesh,
                                       DenseMatrix& Ke, Vector& Fe) const {
    const Element& elem = mesh.element(elem_id);
    ElementType elem_type = elem.type();
    
    // 获取节点坐标
    std::vector<Vec3> coords = get_element_coords(elem_id, mesh);
    Index n_nodes = coords.size();
    
    // 创建形函数对象
    auto shape_func = shape::ShapeFunctionFactory::create(elem_type);
    if (!shape_func) {
        FEM_WARN("ElasticityUnified: unsupported element type");
        return;
    }
    
    // 检查维度一致性
    int elem_dim = shape_func->dimension();
    
    if (dimension_ != elem_dim) {
        FEM_WARN("ElasticityUnified: dimension mismatch (physics: " 
                 + std::to_string(dimension_) + ", element: " 
                 + std::to_string(elem_dim) + ")");
        return;
    }
    
    int dofs_per_node = dimension_;
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = get_gauss_order(elem_type);
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    // 初始化单元矩阵
    Index n_dofs = n_nodes * dofs_per_node;
    Ke.resize(n_dofs, n_dofs);
    Ke.zero();
    Fe.resize(n_dofs);
    for (Index i = 0; i < n_dofs; ++i) {
        Fe[i] = 0.0;
    }
    
    // 创建状态变量（弹性材料为空，塑性材料需要）
    constitutive::StateVariables state = material_->createState();
    
    // 高斯积分循环
    for (size_t gp = 0; gp < gauss_points.size(); ++gp) {
        const Vec3& xi = gauss_points[gp];
        Real w = weights[gp];
        
        // 计算形函数物理坐标导数
        DenseMatrix dN_dxyz;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dxyz);
        
        // 计算雅可比行列式
        DenseMatrix J = shape_func->computeJacobian(xi, coords);
        Real det_J = compute_jacobian_determinant(J);
        
        if (std::abs(det_J) < 1e-15) {
            FEM_WARN("ElasticityUnified: degenerate element (det(J) near zero)");
            continue;
        }
        
        Real dV = w * std::abs(det_J);
        
        // 构造 B 矩阵
        DenseMatrix B;
        if (dimension_ == 2) {
            B = build_B_matrix_2D(dN_dxyz);
        } else {
            B = build_B_matrix_3D(dN_dxyz);
        }
        
        // 获取材料切线刚度矩阵 D
        // 注意：对于线性弹性，D 不依赖应变，但对于非线性材料可能需要当前应变状态
        Vector strain_dummy;  // 弹性材料不需要应变信息
        DenseMatrix D;
        material_->computeTangent(strain_dummy, D, state);
        
        // 装配刚度矩阵: Ke += B^T * D * B * dV
        DenseMatrix DB = D * B;  // D * B (strain_size × n_dofs)
        
        for (Index i = 0; i < n_dofs; ++i) {
            for (Index j = 0; j < n_dofs; ++j) {
                Real sum = 0.0;
                Index strain_size = D.rows();
                for (Index k = 0; k < strain_size; ++k) {
                    sum += B(k, i) * DB(k, j);
                }
                Ke(i, j) += sum * dV;
            }
        }
        
        // Fe = 0 (体力暂不考虑)
    }
}

void ElasticityUnified::compute_stiffness(Index elem_id, const Mesh& mesh,
                                         DenseMatrix& Ke) const {
    Vector Fe_dummy;
    compute_element(elem_id, mesh, Ke, Fe_dummy);
}

// ═══════════════════════════════════════════════════════════
// 访问器实现（仅对简单构造函数有效）
// ═══════════════════════════════════════════════════════════

Real ElasticityUnified::youngs_modulus() const {
    if (!own_material_) {
        FEM_WARN("ElasticityUnified::youngs_modulus() called on custom material");
        return 0.0;
    }
    // 内部创建的一定是 IsotropicElastic
    auto* iso_elastic = dynamic_cast<constitutive::IsotropicElastic*>(material_);
    if (iso_elastic) {
        // IsotropicElastic 有参数访问接口，但我们需要从 elasticityTensor 反推
        // 实际上 IsotropicElastic 应该提供 E() 和 nu() 访问器
        // 这里暂时返回通过拉梅常数反推的结果
        Real mu = iso_elastic->mu();
        Real lambda = iso_elastic->lambda();
        
        // E = mu * (3*lambda + 2*mu) / (lambda + mu)
        Real E = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
        return E;
    }
    return 0.0;
}

Real ElasticityUnified::poissons_ratio() const {
    if (!own_material_) {
        FEM_WARN("ElasticityUnified::poissons_ratio() called on custom material");
        return 0.0;
    }
    auto* iso_elastic = dynamic_cast<constitutive::IsotropicElastic*>(material_);
    if (iso_elastic) {
        Real mu = iso_elastic->mu();
        Real lambda = iso_elastic->lambda();
        
        // nu = lambda / (2 * (lambda + mu))
        Real nu = lambda / (2.0 * (lambda + mu));
        return nu;
    }
    return 0.0;
}

PlaneType ElasticityUnified::plane_type() const {
    if (dimension_ != 2) {
        FEM_WARN("ElasticityUnified::plane_type() called on 3D physics");
    }
    return plane_type_;
}

} // namespace physics
} // namespace fem
