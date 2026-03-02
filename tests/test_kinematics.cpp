#include <gtest/gtest.h>
#include "core/kinematics.h"
#include <cmath>

using namespace fem;

// ═══════════════════════════════════════════════════════════
// 变形梯度测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, DeformationGradient_Identity) {
    // 零位移：∇u = 0 → F = I
    DenseMatrix grad_u(3, 3);
    grad_u.fill(0.0);
    
    DenseMatrix F = Kinematics::deformationGradient(grad_u);
    
    // F 应该是单位矩阵
    EXPECT_NEAR(F(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.0, 1e-14);
    EXPECT_NEAR(F(2, 2), 1.0, 1e-14);
    EXPECT_NEAR(F(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(F(0, 2), 0.0, 1e-14);
    EXPECT_NEAR(F(1, 2), 0.0, 1e-14);
}

TEST(KinematicsTest, DeformationGradient_UniformStretch) {
    // 均匀拉伸：u = 0.1 * x → ∂u/∂x = 0.1 → F_xx = 1.1
    DenseMatrix grad_u(3, 3);
    grad_u.fill(0.0);
    grad_u(0, 0) = 0.1;  // ∂u/∂x
    
    DenseMatrix F = Kinematics::deformationGradient(grad_u);
    
    EXPECT_NEAR(F(0, 0), 1.1, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.0, 1e-14);
    EXPECT_NEAR(F(2, 2), 1.0, 1e-14);
}

TEST(KinematicsTest, DeformationGradient_Shear) {
    // 简单剪切：v = 0.1 * x → ∂v/∂x = 0.1 → F_yx = 0.1
    DenseMatrix grad_u(3, 3);
    grad_u.fill(0.0);
    grad_u(1, 0) = 0.1;  // ∂v/∂x
    
    DenseMatrix F = Kinematics::deformationGradient(grad_u);
    
    EXPECT_NEAR(F(1, 0), 0.1, 1e-14);
    EXPECT_NEAR(F(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.0, 1e-14);
}

TEST(KinematicsTest, DeformationGradientFromVoigt_3D) {
    // 3D Voigt: [∂u/∂x, ∂v/∂y, ∂w/∂z, ∂u/∂y, ∂v/∂x, ∂v/∂z, ∂w/∂y, ∂w/∂x, ∂u/∂z]
    Vector grad_u_voigt(9);
    grad_u_voigt[0] = 0.1;   // ∂u/∂x
    grad_u_voigt[1] = 0.2;   // ∂v/∂y
    grad_u_voigt[2] = 0.3;   // ∂w/∂z
    grad_u_voigt[3] = 0.05;  // ∂u/∂y
    grad_u_voigt[4] = 0.06;  // ∂v/∂x
    grad_u_voigt[5] = 0.07;  // ∂v/∂z
    grad_u_voigt[6] = 0.08;  // ∂w/∂y
    grad_u_voigt[7] = 0.09;  // ∂w/∂x
    grad_u_voigt[8] = 0.04;  // ∂u/∂z
    
    DenseMatrix F = Kinematics::deformationGradientFromVoigt(grad_u_voigt, 3);
    
    EXPECT_NEAR(F(0, 0), 1.1, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.2, 1e-14);
    EXPECT_NEAR(F(2, 2), 1.3, 1e-14);
    EXPECT_NEAR(F(0, 1), 0.05, 1e-14);
    EXPECT_NEAR(F(1, 0), 0.06, 1e-14);
}

TEST(KinematicsTest, DeformationGradientFromVoigt_2D) {
    // 2D Voigt: [∂u/∂x, ∂v/∂y, ∂u/∂y, ∂v/∂x]
    Vector grad_u_voigt(4);
    grad_u_voigt[0] = 0.1;   // ∂u/∂x
    grad_u_voigt[1] = 0.2;   // ∂v/∂y
    grad_u_voigt[2] = 0.05;  // ∂u/∂y
    grad_u_voigt[3] = 0.06;  // ∂v/∂x
    
    DenseMatrix F = Kinematics::deformationGradientFromVoigt(grad_u_voigt, 2);
    
    EXPECT_NEAR(F(0, 0), 1.1, 1e-14);
    EXPECT_NEAR(F(1, 1), 1.2, 1e-14);
    EXPECT_NEAR(F(2, 2), 1.0, 1e-14);  // 平面应变
    EXPECT_NEAR(F(0, 1), 0.05, 1e-14);
    EXPECT_NEAR(F(1, 0), 0.06, 1e-14);
}

// ═══════════════════════════════════════════════════════════
// 应变度量测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, RightCauchyGreen_Identity) {
    // F = I → C = I
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix C = Kinematics::rightCauchyGreen(F);
    
    EXPECT_NEAR(C(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(C(1, 1), 1.0, 1e-14);
    EXPECT_NEAR(C(2, 2), 1.0, 1e-14);
    EXPECT_NEAR(C(0, 1), 0.0, 1e-14);
}

TEST(KinematicsTest, RightCauchyGreen_Stretch) {
    // F = 2I → C = 4I
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 2.0;
    F(2, 2) = 2.0;
    
    DenseMatrix C = Kinematics::rightCauchyGreen(F);
    
    EXPECT_NEAR(C(0, 0), 4.0, 1e-14);
    EXPECT_NEAR(C(1, 1), 4.0, 1e-14);
    EXPECT_NEAR(C(2, 2), 4.0, 1e-14);
}

TEST(KinematicsTest, GreenLagrangeStrain_Identity) {
    // F = I → E = 0
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix E = Kinematics::greenLagrangeStrain(F);
    
    EXPECT_NEAR(E(0, 0), 0.0, 1e-14);
    EXPECT_NEAR(E(1, 1), 0.0, 1e-14);
    EXPECT_NEAR(E(2, 2), 0.0, 1e-14);
    EXPECT_NEAR(E(0, 1), 0.0, 1e-14);
}

TEST(KinematicsTest, GreenLagrangeStrain_SmallStretch) {
    // 小变形：F_xx = 1.01 → E_xx ≈ 0.01（线性近似）
    // 精确：E_xx = 0.5(1.01^2 - 1) = 0.5 * 0.0201 = 0.01005
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.01;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix E = Kinematics::greenLagrangeStrain(F);
    
    EXPECT_NEAR(E(0, 0), 0.01005, 1e-14);
    EXPECT_NEAR(E(1, 1), 0.0, 1e-14);
    EXPECT_NEAR(E(2, 2), 0.0, 1e-14);
}

TEST(KinematicsTest, GreenLagrangeStrain_LargeStretch) {
    // 大变形：F_xx = 2 → E_xx = 0.5(4 - 1) = 1.5
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix E = Kinematics::greenLagrangeStrain(F);
    
    EXPECT_NEAR(E(0, 0), 1.5, 1e-14);
}

TEST(KinematicsTest, StrainToVoigt_3D) {
    DenseMatrix E(3, 3);
    E.fill(0.0);
    E(0, 0) = 0.01;
    E(1, 1) = 0.02;
    E(2, 2) = 0.03;
    E(0, 1) = E(1, 0) = 0.005;
    E(1, 2) = E(2, 1) = 0.006;
    E(0, 2) = E(2, 0) = 0.007;
    
    Vector E_voigt = Kinematics::strainToVoigt(E, 3);
    
    EXPECT_EQ(E_voigt.size(), 6);
    EXPECT_NEAR(E_voigt[0], 0.01, 1e-14);
    EXPECT_NEAR(E_voigt[1], 0.02, 1e-14);
    EXPECT_NEAR(E_voigt[2], 0.03, 1e-14);
    EXPECT_NEAR(E_voigt[3], 0.01, 1e-14);   // 2 * 0.005 = 0.01（工程应变）
    EXPECT_NEAR(E_voigt[4], 0.012, 1e-14);  // 2 * 0.006
    EXPECT_NEAR(E_voigt[5], 0.014, 1e-14);  // 2 * 0.007
}

TEST(KinematicsTest, VoigtToStrain_3D) {
    Vector E_voigt(6);
    E_voigt[0] = 0.01;
    E_voigt[1] = 0.02;
    E_voigt[2] = 0.03;
    E_voigt[3] = 0.01;   // 工程应变 2E_xy
    E_voigt[4] = 0.012;
    E_voigt[5] = 0.014;
    
    DenseMatrix E = Kinematics::voigtToStrain(E_voigt, 3);
    
    EXPECT_NEAR(E(0, 0), 0.01, 1e-14);
    EXPECT_NEAR(E(1, 1), 0.02, 1e-14);
    EXPECT_NEAR(E(2, 2), 0.03, 1e-14);
    EXPECT_NEAR(E(0, 1), 0.005, 1e-14);  // 0.01 / 2
    EXPECT_NEAR(E(1, 0), 0.005, 1e-14);  // 对称
}

// ═══════════════════════════════════════════════════════════
// Jacobian 测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, Jacobian_Identity) {
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    Real J = Kinematics::jacobian(F);
    
    EXPECT_NEAR(J, 1.0, 1e-14);
}

TEST(KinematicsTest, Jacobian_UniformStretch) {
    // 均匀拉伸：F = λI → J = λ^3
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 2.0;
    F(2, 2) = 2.0;
    
    Real J = Kinematics::jacobian(F);
    
    EXPECT_NEAR(J, 8.0, 1e-14);
}

TEST(KinematicsTest, Jacobian_Incompressible) {
    // 近似不可压：F_xx * F_yy * F_zz ≈ 1
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 0.707;  // sqrt(0.5)
    F(2, 2) = 0.707;
    
    Real J = Kinematics::jacobian(F);
    
    EXPECT_NEAR(J, 1.0, 1e-2);  // 约 1.0
}

// ═══════════════════════════════════════════════════════════
// 应力转换测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, PushForwardStress_Identity) {
    // F = I, S → σ = S
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix S(3, 3);
    S.fill(0.0);
    S(0, 0) = 100.0;
    S(1, 1) = 50.0;
    
    DenseMatrix sigma = Kinematics::pushForwardStress(F, S);
    
    EXPECT_NEAR(sigma(0, 0), 100.0, 1e-10);
    EXPECT_NEAR(sigma(1, 1), 50.0, 1e-10);
}

TEST(KinematicsTest, PullBackStress_Identity) {
    // F = I, σ → S = σ
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix sigma(3, 3);
    sigma.fill(0.0);
    sigma(0, 0) = 100.0;
    sigma(1, 1) = 50.0;
    
    DenseMatrix S = Kinematics::pullBackStress(F, sigma);
    
    EXPECT_NEAR(S(0, 0), 100.0, 1e-10);
    EXPECT_NEAR(S(1, 1), 50.0, 1e-10);
}

TEST(KinematicsTest, StressPushPullConsistency) {
    // σ → S → σ 应该保持不变
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.5;
    F(1, 1) = 1.2;
    F(2, 2) = 0.9;
    
    DenseMatrix sigma_original(3, 3);
    sigma_original.fill(0.0);
    sigma_original(0, 0) = 100.0;
    sigma_original(1, 1) = 50.0;
    sigma_original(2, 2) = 75.0;
    
    // σ → S
    DenseMatrix S = Kinematics::pullBackStress(F, sigma_original);
    
    // S → σ
    DenseMatrix sigma_recovered = Kinematics::pushForwardStress(F, S);
    
    // 验证一致性
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(sigma_recovered(i, j), sigma_original(i, j), 1e-10);
        }
    }
}

// ═══════════════════════════════════════════════════════════
// 不变量测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, FirstInvariant) {
    DenseMatrix C(3, 3);
    C.fill(0.0);
    C(0, 0) = 1.5;
    C(1, 1) = 2.0;
    C(2, 2) = 2.5;
    
    Real I1 = Kinematics::firstInvariant(C);
    
    EXPECT_NEAR(I1, 6.0, 1e-14);
}

TEST(KinematicsTest, ThirdInvariant) {
    // I_3 = det(C)，对于 C = F^T F，I_3 = J^2
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 2.0;
    F(2, 2) = 2.0;
    
    DenseMatrix C = Kinematics::rightCauchyGreen(F);
    Real I3 = Kinematics::thirdInvariant(C);
    
    Real J = Kinematics::jacobian(F);
    
    EXPECT_NEAR(I3, J * J, 1e-10);  // I_3 = J^2
}

// ═══════════════════════════════════════════════════════════
// 辅助函数测试
// ═══════════════════════════════════════════════════════════

TEST(KinematicsTest, Det3x3_Identity) {
    DenseMatrix I(3, 3);
    I.fill(0.0);
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    Real det = Kinematics::det3x3(I);
    
    EXPECT_NEAR(det, 1.0, 1e-14);
}

TEST(KinematicsTest, Det3x3_Diagonal) {
    DenseMatrix A(3, 3);
    A.fill(0.0);
    A(0, 0) = 2.0;
    A(1, 1) = 3.0;
    A(2, 2) = 4.0;
    
    Real det = Kinematics::det3x3(A);
    
    EXPECT_NEAR(det, 24.0, 1e-14);  // 2 * 3 * 4
}

TEST(KinematicsTest, Inverse3x3_Identity) {
    DenseMatrix I(3, 3);
    I.fill(0.0);
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    DenseMatrix I_inv = Kinematics::inverse3x3(I);
    
    EXPECT_NEAR(I_inv(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(I_inv(1, 1), 1.0, 1e-14);
    EXPECT_NEAR(I_inv(2, 2), 1.0, 1e-14);
}

TEST(KinematicsTest, Inverse3x3_Diagonal) {
    DenseMatrix A(3, 3);
    A.fill(0.0);
    A(0, 0) = 2.0;
    A(1, 1) = 3.0;
    A(2, 2) = 4.0;
    
    DenseMatrix A_inv = Kinematics::inverse3x3(A);
    
    EXPECT_NEAR(A_inv(0, 0), 0.5, 1e-14);
    EXPECT_NEAR(A_inv(1, 1), 1.0/3.0, 1e-14);
    EXPECT_NEAR(A_inv(2, 2), 0.25, 1e-14);
}

TEST(KinematicsTest, Inverse3x3_VerifyInverse) {
    // A * A_inv = I
    DenseMatrix A(3, 3);
    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0; A(1, 2) = 4.0;
    A(2, 0) = 5.0; A(2, 1) = 6.0; A(2, 2) = 0.0;
    
    DenseMatrix A_inv = Kinematics::inverse3x3(A);
    
    DenseMatrix I = A * A_inv;
    
    // 验证单位矩阵
    EXPECT_NEAR(I(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(I(1, 1), 1.0, 1e-10);
    EXPECT_NEAR(I(2, 2), 1.0, 1e-10);
    EXPECT_NEAR(I(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(I(0, 2), 0.0, 1e-10);
    EXPECT_NEAR(I(1, 2), 0.0, 1e-10);
}
