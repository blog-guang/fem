#include "isotropic_elastic.h"
#include <cmath>
#include <stdexcept>

namespace fem {
namespace constitutive {

IsotropicElastic::IsotropicElastic(Real E, Real nu, 
                                   int dimension, 
                                   bool plane_stress)
    : Material(dimension == 3 ? 6 : 3),
      dimension_(dimension),
      plane_stress_(plane_stress)
{
    setParameter("E", E);
    setParameter("nu", nu);
    validateParameters();
}

void IsotropicElastic::computeStress(
    const Vector& strain_inc, 
    Vector& stress, 
    StateVariables& /*state*/
) {
    // 弹性材料：σ_new = σ_old + D : Δε
    DenseMatrix D = buildElasticityTensor();
    
    Vector delta_stress(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        delta_stress[i] = 0.0;
        for (std::size_t j = 0; j < strain_size_; ++j) {
            delta_stress[i] += D(i, j) * strain_inc[j];
        }
    }
    
    // 更新应力
    if (stress.size() != strain_size_) {
        stress.resize(strain_size_, 0.0);
    }
    stress += delta_stress;
}

void IsotropicElastic::computeTangent(
    const Vector& /*strain*/,
    DenseMatrix& D_mat,
    const StateVariables& /*state*/
) {
    // 弹性材料：切线刚度 = 弹性刚度（常数）
    D_mat = buildElasticityTensor();
}

Real IsotropicElastic::strainEnergy(
    const Vector& strain,
    const StateVariables& /*state*/
) const {
    // Ψ = 0.5 * ε : D : ε
    DenseMatrix D = buildElasticityTensor();
    
    Real energy = 0.0;
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            energy += 0.5 * strain[i] * D(i, j) * strain[j];
        }
    }
    return energy;
}

StateVariables IsotropicElastic::createState() const {
    // 弹性材料：无内部状态变量
    return StateVariables(strain_size_);
}

std::string IsotropicElastic::typeName() const {
    std::string name = "IsotropicElastic";
    if (dimension_ == 2) {
        name += plane_stress_ ? " (plane stress)" : " (plane strain)";
    } else {
        name += " (3D)";
    }
    return name;
}

void IsotropicElastic::validateParameters() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    
    if (E <= 0.0) {
        throw std::invalid_argument("Young's modulus E must be positive");
    }
    
    if (dimension_ == 3) {
        if (nu <= -1.0 || nu >= 0.5) {
            throw std::invalid_argument("Poisson's ratio must be in (-1, 0.5) for 3D");
        }
    } else {
        if (plane_stress_) {
            if (nu <= -1.0 || nu >= 1.0) {
                throw std::invalid_argument("Poisson's ratio must be in (-1, 1) for plane stress");
            }
        } else {
            if (nu <= -1.0 || nu >= 0.5) {
                throw std::invalid_argument("Poisson's ratio must be in (-1, 0.5) for plane strain");
            }
        }
    }
}

Real IsotropicElastic::lambda() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

Real IsotropicElastic::mu() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    return E / (2.0 * (1.0 + nu));
}

Real IsotropicElastic::bulkModulus() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    return E / (3.0 * (1.0 - 2.0 * nu));
}

DenseMatrix IsotropicElastic::elasticityTensor() const {
    return buildElasticityTensor();
}

DenseMatrix IsotropicElastic::buildElasticityTensor() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    
    if (dimension_ == 3) {
        // 3D弹性刚度矩阵（6×6，Voigt记号）
        Real lambda_val = lambda();
        Real mu_val = mu();
        
        DenseMatrix D(6, 6);
        D.zero();
        
        // 对角块（正应力）
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                D(i, j) = lambda_val;
            }
            D(i, i) += 2.0 * mu_val;
        }
        
        // 剪应力
        D(3, 3) = mu_val;  // σ12
        D(4, 4) = mu_val;  // σ23
        D(5, 5) = mu_val;  // σ13
        
        return D;
        
    } else {
        // 2D弹性刚度矩阵（3×3）
        DenseMatrix D(3, 3);
        D.zero();
        
        if (plane_stress_) {
            // 平面应力
            Real factor = E / (1.0 - nu * nu);
            D(0, 0) = factor;
            D(0, 1) = factor * nu;
            D(1, 0) = factor * nu;
            D(1, 1) = factor;
            D(2, 2) = factor * (1.0 - nu) / 2.0;
        } else {
            // 平面应变
            Real factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            D(0, 0) = factor * (1.0 - nu);
            D(0, 1) = factor * nu;
            D(1, 0) = factor * nu;
            D(1, 1) = factor * (1.0 - nu);
            D(2, 2) = factor * (1.0 - 2.0 * nu) / 2.0;
        }
        
        return D;
    }
}

}  // namespace constitutive
}  // namespace fem
