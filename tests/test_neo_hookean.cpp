#include <gtest/gtest.h>
#include "material/neo_hookean.h"
#include "material/material_factory.h"
#include <cmath>

using namespace fem;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════
// NeoHookean 基础测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, Construction_C10D1) {
    Real C10 = 1e6;  // 1 MPa
    Real D1 = 1e-6;
    
    NeoHookean mat(C10, D1, 3);
    
    EXPECT_EQ(mat.typeName(), "NeoHookean_3D");
}

TEST(NeoHookeanTest, Construction_Engineering) {
    Real E = 10e6;   // 10 MPa (橡胶级别)
    Real nu = 0.45;  // 近似不可压
    
    NeoHookean mat(E, nu, 3, true);
    
    EXPECT_EQ(mat.typeName(), "NeoHookean_3D");
}

TEST(NeoHookeanTest, InvalidParameters) {
    // C10 <= 0
    EXPECT_THROW(NeoHookean(-1.0, 1e-6, 3), std::invalid_argument);
    
    // D1 <= 0
    EXPECT_THROW(NeoHookean(1e6, -1e-6, 3), std::invalid_argument);
    
    // nu >= 0.5 (不可压)
    EXPECT_THROW(NeoHookean(10e6, 0.5, 3, true), std::invalid_argument);
    
    // 错误维度
    EXPECT_THROW(NeoHookean(1e6, 1e-6, 1), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════
// 小变形测试（线弹性近似）
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, SmallStrain_3D_Uniaxial) {
    // 等效材料参数
    Real E = 10e6;   // 10 MPa
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    // 单轴拉伸应变：ε_xx = 0.001（0.1%应变，小变形）
    Vector strain(6, 0.0);
    strain[0] = 0.001;  // ε_xx
    
    Vector stress;
    mat.computeStress(strain, stress, state);
    
    // 验证应力大小（线弹性近似）
    EXPECT_GT(stress[0], 0.0);  // 拉应力为正
    EXPECT_GT(stress[0], 5000.0);  // 至少 5 kPa
    EXPECT_LT(stress[0], 20000.0); // 小于 20 kPa
    
    // 注意：在约束条件下（ε_yy = ε_zz = 0），会有横向应力
    // σ_yy = σ_zz = ν/(1-ν) * σ_xx（泊松效应）
    // 这里测试 σ_yy, σ_zz 相等
    EXPECT_NEAR(stress[1], stress[2], 1e-6);
    
    // 剪切应力应该为 0
    EXPECT_NEAR(stress[3], 0.0, 1e-6);
    EXPECT_NEAR(stress[4], 0.0, 1e-6);
    EXPECT_NEAR(stress[5], 0.0, 1e-6);
}

TEST(NeoHookeanTest, SmallStrain_3D_PureShear) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    // 纯剪切：γ_xy = 0.001
    Vector strain(6, 0.0);
    strain[3] = 0.001;  // γ_xy
    
    Vector stress;
    mat.computeStress(strain, stress, state);
    
    // 剪应力：τ_xy = G * γ_xy，其中 G = E / (2(1+ν))
    Real G = E / (2.0 * (1.0 + nu));
    Real expected_shear = G * strain[3];
    
    EXPECT_NEAR(stress[3], expected_shear, 1e-6);
    
    // 法向应力应该为 0
    EXPECT_NEAR(stress[0], 0.0, 1e-6);
    EXPECT_NEAR(stress[1], 0.0, 1e-6);
    EXPECT_NEAR(stress[2], 0.0, 1e-6);
}

TEST(NeoHookeanTest, SmallStrain_2D_PlaneStrain) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 2, true);
    StateVariables state = mat.createState();
    
    // 2D 平面应变单轴拉伸
    Vector strain(4, 0.0);
    strain[0] = 0.001;  // ε_xx
    
    Vector stress;
    mat.computeStress(strain, stress, state);
    
    EXPECT_GT(stress[0], 0.0);
    
    // 平面应变：σ_zz ≠ 0
    EXPECT_NE(stress[2], 0.0);
}

// ═══════════════════════════════════════════════════════════
// 切线刚度矩阵测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, TangentStiffness_3D) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    Vector strain(6, 0.0);
    DenseMatrix D;
    
    mat.computeTangent(strain, D, state);
    
    // 检查维度
    EXPECT_EQ(D.rows(), 6);
    EXPECT_EQ(D.cols(), 6);
    
    // 检查对称性
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            EXPECT_NEAR(D(i, j), D(j, i), 1e-10);
        }
    }
    
    // 检查正定性（对角元素为正）
    for (int i = 0; i < 6; i++) {
        EXPECT_GT(D(i, i), 0.0);
    }
}

TEST(NeoHookeanTest, TangentStiffness_2D) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 2, true);
    StateVariables state = mat.createState();
    
    Vector strain(4, 0.0);
    DenseMatrix D;
    
    mat.computeTangent(strain, D, state);
    
    EXPECT_EQ(D.rows(), 4);
    EXPECT_EQ(D.cols(), 4);
    
    // 对称性
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_NEAR(D(i, j), D(j, i), 1e-10);
        }
    }
}

// ═══════════════════════════════════════════════════════════
// 应变能测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, StrainEnergy) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    // 单轴拉伸
    Vector strain(6, 0.0);
    strain[0] = 0.001;
    
    Real W = mat.strainEnergy(strain, state);
    
    // 应变能应该为正
    EXPECT_GT(W, 0.0);
    
    // 应变能应该随应变增加而增加
    strain[0] = 0.002;
    Real W2 = mat.strainEnergy(strain, state);
    EXPECT_GT(W2, W);
}

TEST(NeoHookeanTest, StrainEnergy_ZeroStrain) {
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    Vector strain(6, 0.0);
    Real W = mat.strainEnergy(strain, state);
    
    EXPECT_NEAR(W, 0.0, 1e-10);
}

// ═══════════════════════════════════════════════════════════
// MaterialFactory 集成测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, Factory_C10D1) {
    auto mat = MaterialFactory::create("NeoHookean", {
        {"C10", 1e6},
        {"D1", 1e-6},
        {"dimension", 3}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_TRUE(mat->typeName().find("NeoHookean") != std::string::npos);
}

TEST(NeoHookeanTest, Factory_Engineering) {
    auto mat = MaterialFactory::create("NeoHookean", {
        {"E", 10e6},
        {"nu", 0.3},
        {"dimension", 3}
    });
    
    ASSERT_NE(mat, nullptr);
    
    // 验证能够计算应力
    StateVariables state = mat->createState();
    Vector strain(6, 0.0);
    strain[0] = 0.001;
    
    Vector stress;
    mat->computeStress(strain, stress, state);
    
    EXPECT_GT(stress[0], 0.0);
}

TEST(NeoHookeanTest, Factory_2D) {
    auto mat = MaterialFactory::create("NeoHookean", {
        {"E", 10e6},
        {"nu", 0.3},
        {"dimension", 2}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_TRUE(mat->typeName().find("2D") != std::string::npos);
}

// ═══════════════════════════════════════════════════════════
// 辅助函数测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, RightCauchyGreen) {
    // 单位变形梯度 F = I
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    DenseMatrix C = NeoHookean::computeRightCauchyGreen(F);
    
    // C = F^T F = I
    EXPECT_NEAR(C(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(C(1, 1), 1.0, 1e-10);
    EXPECT_NEAR(C(2, 2), 1.0, 1e-10);
    EXPECT_NEAR(C(0, 1), 0.0, 1e-10);
}

TEST(NeoHookeanTest, FirstInvariant) {
    DenseMatrix C(3, 3);
    C.fill(0.0);
    C(0, 0) = 1.5;
    C(1, 1) = 2.0;
    C(2, 2) = 2.5;
    
    Real I1 = NeoHookean::firstInvariant(C);
    
    EXPECT_NEAR(I1, 6.0, 1e-10);  // 1.5 + 2.0 + 2.5
}

TEST(NeoHookeanTest, Jacobian_Identity) {
    // F = I → J = 1
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    Real J = NeoHookean::jacobian(F);
    
    EXPECT_NEAR(J, 1.0, 1e-10);
}

TEST(NeoHookeanTest, Jacobian_Scaled) {
    // 均匀拉伸：F = 2I → J = 8
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 2.0;
    F(1, 1) = 2.0;
    F(2, 2) = 2.0;
    
    Real J = NeoHookean::jacobian(F);
    
    EXPECT_NEAR(J, 8.0, 1e-10);
}

// ═══════════════════════════════════════════════════════════
// 大变形测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, LargeDeformation_2ndPiolaKirchhoff) {
    // 测试 2nd Piola-Kirchhoff 应力计算
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    
    // 单位变形梯度 F = I → S = 0（无应力）
    DenseMatrix F_identity(3, 3);
    F_identity.fill(0.0);
    F_identity(0, 0) = 1.0;
    F_identity(1, 1) = 1.0;
    F_identity(2, 2) = 1.0;
    
    DenseMatrix S = mat.compute2ndPiolaKirchhoff(F_identity);
    
    // 单位变形时应力应该接近零
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(S(i, j), 0.0, 1e-6);
        }
    }
}

TEST(NeoHookeanTest, LargeDeformation_UniformStretch) {
    // 测试均匀拉伸
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    
    // 均匀拉伸 10%: F_ii = 1.1
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.1;
    F(1, 1) = 1.1;
    F(2, 2) = 1.1;
    
    DenseMatrix S = mat.compute2ndPiolaKirchhoff(F);
    
    // 均匀变形时，主应力应该相等
    EXPECT_NEAR(S(0, 0), S(1, 1), 1e-6);
    EXPECT_NEAR(S(1, 1), S(2, 2), 1e-6);
    
    // 剪切应力应该为零
    EXPECT_NEAR(S(0, 1), 0.0, 1e-8);
    EXPECT_NEAR(S(0, 2), 0.0, 1e-8);
    EXPECT_NEAR(S(1, 2), 0.0, 1e-8);
}

TEST(NeoHookeanTest, LargeDeformation_Uniaxial) {
    // 单轴拉伸测试
    
    Real E = 10e6;
    Real nu = 0.49;  // 近不可压
    
    NeoHookean mat(E, nu, 3, true);
    
    // 单轴拉伸 20%: F_xx = 1.2
    // 不可压条件：F_yy * F_zz ≈ 1/F_xx ≈ 0.833
    Real lambda_x = 1.2;
    Real lambda_yz = 1.0 / std::sqrt(lambda_x);  // ≈ 0.913
    
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = lambda_x;
    F(1, 1) = lambda_yz;
    F(2, 2) = lambda_yz;
    
    DenseMatrix S = mat.compute2ndPiolaKirchhoff(F);
    
    // 拉伸方向应力应该为正
    EXPECT_GT(S(0, 0), 0.0);
    
    // 横向应力应该接近（考虑泊松效应）
    EXPECT_NEAR(S(1, 1), S(2, 2), 1e-6);
}

TEST(NeoHookeanTest, LargeDeformation_CauchyStress) {
    // 测试 Cauchy 应力计算
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    
    // 单位变形
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    Vector sigma_voigt;
    mat.computeStressFromF(F, sigma_voigt);
    
    EXPECT_EQ(sigma_voigt.size(), 6);
    
    // 单位变形时 Cauchy 应力应该接近零
    for (int i = 0; i < 6; i++) {
        EXPECT_NEAR(sigma_voigt[i], 0.0, 1e-6);
    }
}

TEST(NeoHookeanTest, LargeDeformation_StressVoigtConversion) {
    // 测试应力张量到 Voigt 记号的转换
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    
    // 拉伸 15%
    DenseMatrix F(3, 3);
    F.fill(0.0);
    F(0, 0) = 1.15;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    
    Vector sigma_voigt;
    mat.computeStressFromF(F, sigma_voigt);
    
    // 拉伸方向应力为正
    EXPECT_GT(sigma_voigt[0], 0.0);
    
    // 剪切应力应该为零（单轴变形）
    EXPECT_NEAR(sigma_voigt[3], 0.0, 1e-6);  // σ_xy
    EXPECT_NEAR(sigma_voigt[4], 0.0, 1e-6);  // σ_yz
    EXPECT_NEAR(sigma_voigt[5], 0.0, 1e-6);  // σ_xz
}

// ═══════════════════════════════════════════════════════════
// 小变形自动识别测试
// ═══════════════════════════════════════════════════════════

TEST(NeoHookeanTest, AutoDetect_SmallStrain) {
    // 小应变应自动使用线弹性近似
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    // 小应变（0.5%）
    Vector strain(6, 0.0);
    strain[0] = 0.005;
    
    Vector stress;
    mat.computeStress(strain, stress, state);
    
    // 应该使用线弹性近似
    EXPECT_GT(stress[0], 0.0);
    
    // 验证线弹性关系（近似）
    Real G = E / (2.0 * (1.0 + nu));
    Real K = E / (3.0 * (1.0 - 2.0 * nu));
    Real lambda = K - 2.0 * G / 3.0;
    
    Real expected_stress = lambda * strain[0] + 2.0 * G * strain[0];
    EXPECT_NEAR(stress[0], expected_stress, 1e-3);
}

TEST(NeoHookeanTest, AutoDetect_LargeStrain) {
    // 大应变应触发大变形公式（当前仍使用线弹性近似）
    
    Real E = 10e6;
    Real nu = 0.3;
    
    NeoHookean mat(E, nu, 3, true);
    StateVariables state = mat.createState();
    
    // 大应变（15%）
    Vector strain(6, 0.0);
    strain[0] = 0.15;
    
    Vector stress;
    mat.computeStress(strain, stress, state);
    
    // 应该有应力响应
    EXPECT_GT(stress[0], 0.0);
    
    // 注意：当前大变形分支仍使用线弹性近似
    // 完整实现需要存储变形梯度 F
}
