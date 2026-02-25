#include "orthotropic_elastic.h"
#include <cmath>
#include <stdexcept>

namespace fem {
namespace constitutive {

OrthotropicElastic::OrthotropicElastic(Real E1, Real E2, Real nu12, Real G12, 
                                       Real theta)
    : Material(3),  // 2D: 3分量
      dimension_(2),
      theta_(theta * M_PI / 180.0)  // 转换为弧度
{
    setParameter("E1", E1);
    setParameter("E2", E2);
    setParameter("nu12", nu12);
    setParameter("G12", G12);
    setParameter("theta", theta);
    validateParameters();
}

OrthotropicElastic::OrthotropicElastic(Real E1, Real E2, Real E3,
                                       Real nu12, Real nu23, Real nu13,
                                       Real G12, Real G23, Real G13)
    : Material(6),  // 3D: 6分量
      dimension_(3),
      theta_(0.0)
{
    setParameter("E1", E1);
    setParameter("E2", E2);
    setParameter("E3", E3);
    setParameter("nu12", nu12);
    setParameter("nu23", nu23);
    setParameter("nu13", nu13);
    setParameter("G12", G12);
    setParameter("G23", G23);
    setParameter("G13", G13);
    validateParameters();
}

void OrthotropicElastic::computeStress(
    const Vector& strain_inc, 
    Vector& stress, 
    StateVariables& /*state*/
) {
    // 弹性材料：σ_new = σ_old + D : Δε
    DenseMatrix D = globalStiffness();
    
    Vector delta_stress(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        delta_stress[i] = 0.0;
        for (std::size_t j = 0; j < strain_size_; ++j) {
            delta_stress[i] += D(i, j) * strain_inc[j];
        }
    }
    
    if (stress.size() != strain_size_) {
        stress.resize(strain_size_, 0.0);
    }
    stress += delta_stress;
}

void OrthotropicElastic::computeTangent(
    const Vector& /*strain*/,
    DenseMatrix& D_mat,
    const StateVariables& /*state*/
) {
    // 弹性材料：切线刚度 = 弹性刚度（常数）
    D_mat = globalStiffness();
}

Real OrthotropicElastic::strainEnergy(
    const Vector& strain,
    const StateVariables& /*state*/
) const {
    // Ψ = 0.5 * ε : D : ε
    DenseMatrix D = globalStiffness();
    
    Real energy = 0.0;
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            energy += 0.5 * strain[i] * D(i, j) * strain[j];
        }
    }
    return energy;
}

StateVariables OrthotropicElastic::createState() const {
    return StateVariables(strain_size_);
}

std::string OrthotropicElastic::typeName() const {
    std::string name = "OrthotropicElastic";
    if (dimension_ == 2) {
        name += " (2D, theta=" + std::to_string(theta_ * 180.0 / M_PI) + "°)";
    } else {
        name += " (3D)";
    }
    return name;
}

void OrthotropicElastic::validateParameters() const {
    Real E1 = getParameter("E1");
    Real E2 = getParameter("E2");
    Real nu12 = getParameter("nu12");
    Real G12 = getParameter("G12");
    
    if (E1 <= 0.0 || E2 <= 0.0) {
        throw std::invalid_argument("Young's moduli E1, E2 must be positive");
    }
    if (G12 <= 0.0) {
        throw std::invalid_argument("Shear modulus G12 must be positive");
    }
    
    // 泊松比约束：|nu12| < sqrt(E1/E2)
    Real nu_max = std::sqrt(E1 / E2);
    if (std::abs(nu12) >= nu_max) {
        throw std::invalid_argument("Poisson's ratio must satisfy |nu12| < sqrt(E1/E2)");
    }
    
    if (dimension_ == 3) {
        Real E3 = getParameter("E3");
        Real nu23 = getParameter("nu23");
        Real nu13 = getParameter("nu13");
        Real G23 = getParameter("G23");
        Real G13 = getParameter("G13");
        
        if (E3 <= 0.0 || G23 <= 0.0 || G13 <= 0.0) {
            throw std::invalid_argument("All moduli must be positive");
        }
        
        // 3D泊松比约束（简化检查）
        if (std::abs(nu23) >= 0.5 || std::abs(nu13) >= 0.5) {
            throw std::invalid_argument("Poisson's ratios must be < 0.5");
        }
    }
}

DenseMatrix OrthotropicElastic::localStiffness() const {
    if (dimension_ == 2) {
        return buildLocalStiffness2D();
    } else {
        return buildLocalStiffness3D();
    }
}

DenseMatrix OrthotropicElastic::globalStiffness() const {
    DenseMatrix D_local = localStiffness();
    
    if (dimension_ == 2 && std::abs(theta_) > 1e-10) {
        // 需要坐标变换
        DenseMatrix T = transformationMatrix2D();
        DenseMatrix T_inv = T.transpose();  // T是正交矩阵
        
        // D_global = T^{-T} D_local T^{-1}
        DenseMatrix temp = D_local * T_inv;
        return T_inv.transpose() * temp;
    }
    
    return D_local;
}

void OrthotropicElastic::setRotation(Real theta_deg) {
    theta_ = theta_deg * M_PI / 180.0;
    setParameter("theta", theta_deg);
}

DenseMatrix OrthotropicElastic::transformationMatrix2D() const {
    Real c = std::cos(theta_);
    Real s = std::sin(theta_);
    Real c2 = c * c;
    Real s2 = s * s;
    Real cs = c * s;
    
    DenseMatrix T(3, 3);
    T.zero();
    
    // 应力变换矩阵（Voigt记号）
    T(0, 0) = c2;
    T(0, 1) = s2;
    T(0, 2) = 2.0 * cs;
    
    T(1, 0) = s2;
    T(1, 1) = c2;
    T(1, 2) = -2.0 * cs;
    
    T(2, 0) = -cs;
    T(2, 1) = cs;
    T(2, 2) = c2 - s2;
    
    return T;
}

DenseMatrix OrthotropicElastic::buildLocalStiffness2D() const {
    Real E1 = getParameter("E1");
    Real E2 = getParameter("E2");
    Real nu12 = getParameter("nu12");
    Real G12 = getParameter("G12");
    
    // nu21 = nu12 * E2 / E1 (互易关系)
    Real nu21 = nu12 * E2 / E1;
    
    // 分母
    Real denom = 1.0 - nu12 * nu21;
    
    DenseMatrix D(3, 3);
    D.zero();
    
    // 平面应力刚度矩阵（主轴坐标系）
    D(0, 0) = E1 / denom;
    D(0, 1) = nu12 * E2 / denom;
    D(1, 0) = nu21 * E1 / denom;
    D(1, 1) = E2 / denom;
    D(2, 2) = G12;
    
    return D;
}

DenseMatrix OrthotropicElastic::buildLocalStiffness3D() const {
    Real E1 = getParameter("E1");
    Real E2 = getParameter("E2");
    Real E3 = getParameter("E3");
    Real nu12 = getParameter("nu12");
    Real nu23 = getParameter("nu23");
    Real nu13 = getParameter("nu13");
    Real G12 = getParameter("G12");
    Real G23 = getParameter("G23");
    Real G13 = getParameter("G13");
    
    // 互易泊松比
    Real nu21 = nu12 * E2 / E1;
    Real nu32 = nu23 * E3 / E2;
    Real nu31 = nu13 * E3 / E1;
    
    // 柔度矩阵元素
    Real gamma = 1.0 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2.0*nu12*nu23*nu31;
    
    DenseMatrix D(6, 6);
    D.zero();
    
    // 主对角块（正应力）
    D(0, 0) = E1 * (1.0 - nu23*nu32) / gamma;
    D(0, 1) = E1 * (nu21 + nu31*nu23) / gamma;
    D(0, 2) = E1 * (nu31 + nu21*nu32) / gamma;
    
    D(1, 0) = E2 * (nu12 + nu32*nu13) / gamma;
    D(1, 1) = E2 * (1.0 - nu13*nu31) / gamma;
    D(1, 2) = E2 * (nu32 + nu12*nu31) / gamma;
    
    D(2, 0) = E3 * (nu13 + nu12*nu23) / gamma;
    D(2, 1) = E3 * (nu23 + nu13*nu21) / gamma;
    D(2, 2) = E3 * (1.0 - nu12*nu21) / gamma;
    
    // 剪应力
    D(3, 3) = G12;
    D(4, 4) = G23;
    D(5, 5) = G13;
    
    return D;
}

}  // namespace constitutive
}  // namespace fem
