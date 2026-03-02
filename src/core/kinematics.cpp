/**
 * kinematics.cpp - 几何非线性运动学实现
 */

#include "core/kinematics.h"
#include <cmath>
#include <stdexcept>

namespace fem {

// ═══════════════════════════════════════════════════════════
// 变形梯度
// ═══════════════════════════════════════════════════════════

DenseMatrix Kinematics::deformationGradient(const DenseMatrix& grad_u) {
    if (grad_u.rows() != 3 || grad_u.cols() != 3) {
        throw std::invalid_argument("grad_u must be 3x3");
    }
    
    // F = I + ∇u
    DenseMatrix F(3, 3);
    F(0, 0) = 1.0 + grad_u(0, 0);
    F(0, 1) = grad_u(0, 1);
    F(0, 2) = grad_u(0, 2);
    
    F(1, 0) = grad_u(1, 0);
    F(1, 1) = 1.0 + grad_u(1, 1);
    F(1, 2) = grad_u(1, 2);
    
    F(2, 0) = grad_u(2, 0);
    F(2, 1) = grad_u(2, 1);
    F(2, 2) = 1.0 + grad_u(2, 2);
    
    return F;
}

DenseMatrix Kinematics::deformationGradientFromVoigt(const Vector& grad_u_voigt, int dimension) {
    DenseMatrix grad_u(3, 3);
    grad_u.fill(0.0);
    
    if (dimension == 3) {
        // 3D: [∂u/∂x, ∂v/∂y, ∂w/∂z, ∂u/∂y, ∂v/∂x, ∂v/∂z, ∂w/∂y, ∂w/∂x, ∂u/∂z]
        if (grad_u_voigt.size() != 9) {
            throw std::invalid_argument("3D grad_u_voigt must have 9 components");
        }
        
        grad_u(0, 0) = grad_u_voigt[0];  // ∂u/∂x
        grad_u(1, 1) = grad_u_voigt[1];  // ∂v/∂y
        grad_u(2, 2) = grad_u_voigt[2];  // ∂w/∂z
        grad_u(0, 1) = grad_u_voigt[3];  // ∂u/∂y
        grad_u(1, 0) = grad_u_voigt[4];  // ∂v/∂x
        grad_u(1, 2) = grad_u_voigt[5];  // ∂v/∂z
        grad_u(2, 1) = grad_u_voigt[6];  // ∂w/∂y
        grad_u(2, 0) = grad_u_voigt[7];  // ∂w/∂x
        grad_u(0, 2) = grad_u_voigt[8];  // ∂u/∂z
    } else {
        // 2D: [∂u/∂x, ∂v/∂y, ∂u/∂y, ∂v/∂x]
        if (grad_u_voigt.size() != 4) {
            throw std::invalid_argument("2D grad_u_voigt must have 4 components");
        }
        
        grad_u(0, 0) = grad_u_voigt[0];  // ∂u/∂x
        grad_u(1, 1) = grad_u_voigt[1];  // ∂v/∂y
        grad_u(0, 1) = grad_u_voigt[2];  // ∂u/∂y
        grad_u(1, 0) = grad_u_voigt[3];  // ∂v/∂x
        // grad_u(2, 2) = 0.0  (平面应变：∂w/∂z = 0)
    }
    
    return deformationGradient(grad_u);
}

// ═══════════════════════════════════════════════════════════
// 应变度量
// ═══════════════════════════════════════════════════════════

DenseMatrix Kinematics::rightCauchyGreen(const DenseMatrix& F) {
    // C = F^T F
    return F.transpose() * F;
}

DenseMatrix Kinematics::greenLagrangeStrain(const DenseMatrix& F) {
    // E = 1/2(F^T F - I)
    DenseMatrix C = rightCauchyGreen(F);
    
    // E = 1/2(C - I)
    DenseMatrix E(3, 3);
    E(0, 0) = 0.5 * (C(0, 0) - 1.0);
    E(0, 1) = 0.5 * C(0, 1);
    E(0, 2) = 0.5 * C(0, 2);
    
    E(1, 0) = 0.5 * C(1, 0);
    E(1, 1) = 0.5 * (C(1, 1) - 1.0);
    E(1, 2) = 0.5 * C(1, 2);
    
    E(2, 0) = 0.5 * C(2, 0);
    E(2, 1) = 0.5 * C(2, 1);
    E(2, 2) = 0.5 * (C(2, 2) - 1.0);
    
    return E;
}

Vector Kinematics::strainToVoigt(const DenseMatrix& E, int dimension) {
    Vector E_voigt;
    
    if (dimension == 3) {
        E_voigt.resize(6);
        E_voigt[0] = E(0, 0);        // E_xx
        E_voigt[1] = E(1, 1);        // E_yy
        E_voigt[2] = E(2, 2);        // E_zz
        E_voigt[3] = 2.0 * E(0, 1);  // 2E_xy（工程应变）
        E_voigt[4] = 2.0 * E(1, 2);  // 2E_yz
        E_voigt[5] = 2.0 * E(0, 2);  // 2E_xz
    } else {
        E_voigt.resize(4);
        E_voigt[0] = E(0, 0);        // E_xx
        E_voigt[1] = E(1, 1);        // E_yy
        E_voigt[2] = E(2, 2);        // E_zz（平面应变）
        E_voigt[3] = 2.0 * E(0, 1);  // 2E_xy
    }
    
    return E_voigt;
}

DenseMatrix Kinematics::voigtToStrain(const Vector& E_voigt, int dimension) {
    DenseMatrix E(3, 3);
    E.fill(0.0);
    
    if (dimension == 3) {
        E(0, 0) = E_voigt[0];
        E(1, 1) = E_voigt[1];
        E(2, 2) = E_voigt[2];
        E(0, 1) = E(1, 0) = E_voigt[3] / 2.0;
        E(1, 2) = E(2, 1) = E_voigt[4] / 2.0;
        E(0, 2) = E(2, 0) = E_voigt[5] / 2.0;
    } else {
        E(0, 0) = E_voigt[0];
        E(1, 1) = E_voigt[1];
        E(2, 2) = E_voigt[2];
        E(0, 1) = E(1, 0) = E_voigt[3] / 2.0;
    }
    
    return E;
}

// ═══════════════════════════════════════════════════════════
// Jacobian
// ═══════════════════════════════════════════════════════════

Real Kinematics::jacobian(const DenseMatrix& F) {
    return det3x3(F);
}

// ═══════════════════════════════════════════════════════════
// 应力转换
// ═══════════════════════════════════════════════════════════

DenseMatrix Kinematics::pushForwardStress(const DenseMatrix& F, const DenseMatrix& S) {
    // σ = J^(-1) F S F^T
    Real J = jacobian(F);
    
    if (J <= 0.0) {
        throw std::runtime_error("Invalid Jacobian (J <= 0): element inversion detected");
    }
    
    DenseMatrix FSF = F * S * F.transpose();
    
    return FSF * (1.0 / J);
}

DenseMatrix Kinematics::pullBackStress(const DenseMatrix& F, const DenseMatrix& sigma) {
    // S = J F^(-1) σ F^(-T)
    Real J = jacobian(F);
    
    if (J <= 0.0) {
        throw std::runtime_error("Invalid Jacobian (J <= 0): element inversion detected");
    }
    
    DenseMatrix F_inv = inverse3x3(F);
    DenseMatrix F_inv_T = F_inv.transpose();
    
    DenseMatrix S = F_inv * sigma * F_inv_T;
    
    return S * J;
}

Vector Kinematics::stressToVoigt(const DenseMatrix& stress, int dimension) {
    Vector stress_voigt;
    
    if (dimension == 3) {
        stress_voigt.resize(6);
        stress_voigt[0] = stress(0, 0);  // σ_xx
        stress_voigt[1] = stress(1, 1);  // σ_yy
        stress_voigt[2] = stress(2, 2);  // σ_zz
        stress_voigt[3] = stress(0, 1);  // σ_xy
        stress_voigt[4] = stress(1, 2);  // σ_yz
        stress_voigt[5] = stress(0, 2);  // σ_xz
    } else {
        stress_voigt.resize(4);
        stress_voigt[0] = stress(0, 0);  // σ_xx
        stress_voigt[1] = stress(1, 1);  // σ_yy
        stress_voigt[2] = stress(2, 2);  // σ_zz
        stress_voigt[3] = stress(0, 1);  // σ_xy
    }
    
    return stress_voigt;
}

DenseMatrix Kinematics::voigtToStress(const Vector& stress_voigt, int dimension) {
    DenseMatrix stress(3, 3);
    stress.fill(0.0);
    
    if (dimension == 3) {
        stress(0, 0) = stress_voigt[0];
        stress(1, 1) = stress_voigt[1];
        stress(2, 2) = stress_voigt[2];
        stress(0, 1) = stress(1, 0) = stress_voigt[3];
        stress(1, 2) = stress(2, 1) = stress_voigt[4];
        stress(0, 2) = stress(2, 0) = stress_voigt[5];
    } else {
        stress(0, 0) = stress_voigt[0];
        stress(1, 1) = stress_voigt[1];
        stress(2, 2) = stress_voigt[2];
        stress(0, 1) = stress(1, 0) = stress_voigt[3];
    }
    
    return stress;
}

// ═══════════════════════════════════════════════════════════
// 不变量
// ═══════════════════════════════════════════════════════════

Real Kinematics::firstInvariant(const DenseMatrix& C) {
    // I_1 = tr(C)
    return C(0, 0) + C(1, 1) + C(2, 2);
}

Real Kinematics::secondInvariant(const DenseMatrix& C) {
    // I_2 = 1/2[(tr C)^2 - tr(C^2)]
    Real tr_C = firstInvariant(C);
    
    // C^2 = C * C
    DenseMatrix C2 = C * C;
    Real tr_C2 = C2(0, 0) + C2(1, 1) + C2(2, 2);
    
    return 0.5 * (tr_C * tr_C - tr_C2);
}

Real Kinematics::thirdInvariant(const DenseMatrix& C) {
    // I_3 = det(C)
    return det3x3(C);
}

// ═══════════════════════════════════════════════════════════
// 辅助函数
// ═══════════════════════════════════════════════════════════

Real Kinematics::det3x3(const DenseMatrix& A) {
    if (A.rows() != 3 || A.cols() != 3) {
        throw std::invalid_argument("Matrix must be 3x3");
    }
    
    // Sarrus 法则
    Real det = A(0,0) * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
             - A(0,1) * (A(1,0)*A(2,2) - A(1,2)*A(2,0))
             + A(0,2) * (A(1,0)*A(2,1) - A(1,1)*A(2,0));
    
    return det;
}

DenseMatrix Kinematics::inverse3x3(const DenseMatrix& A) {
    if (A.rows() != 3 || A.cols() != 3) {
        throw std::invalid_argument("Matrix must be 3x3");
    }
    
    Real det = det3x3(A);
    
    if (std::abs(det) < 1e-14) {
        throw std::runtime_error("Matrix is singular (det ≈ 0)");
    }
    
    // 伴随矩阵法求逆：A^(-1) = adj(A) / det(A)
    // adj(A) = 代数余子式矩阵的转置
    
    DenseMatrix cofactor(3, 3);
    
    // 计算代数余子式矩阵（带符号的余子式）
    // cofactor(i, j) = (-1)^(i+j) * minor(i, j)
    cofactor(0, 0) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1));
    cofactor(0, 1) = -(A(1,0)*A(2,2) - A(1,2)*A(2,0));
    cofactor(0, 2) =  (A(1,0)*A(2,1) - A(1,1)*A(2,0));
    
    cofactor(1, 0) = -(A(0,1)*A(2,2) - A(0,2)*A(2,1));
    cofactor(1, 1) =  (A(0,0)*A(2,2) - A(0,2)*A(2,0));
    cofactor(1, 2) = -(A(0,0)*A(2,1) - A(0,1)*A(2,0));
    
    cofactor(2, 0) =  (A(0,1)*A(1,2) - A(0,2)*A(1,1));
    cofactor(2, 1) = -(A(0,0)*A(1,2) - A(0,2)*A(1,0));
    cofactor(2, 2) =  (A(0,0)*A(1,1) - A(0,1)*A(1,0));
    
    // 伴随矩阵 = 代数余子式矩阵的转置
    DenseMatrix adj = cofactor.transpose();
    
    // A^(-1) = adj(A) / det(A)
    Real inv_det = 1.0 / det;
    DenseMatrix A_inv = adj * inv_det;
    
    return A_inv;
}

}  // namespace fem
