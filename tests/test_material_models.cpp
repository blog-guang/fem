#include <gtest/gtest.h>
#include "material/isotropic_elastic.h"
#include "material/j2_plasticity.h"
#include <cmath>

using namespace fem;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════════
// IsotropicElastic 测试
// ═══════════════════════════════════════════════════════════════

TEST(IsotropicElasticTest, Construction) {
    IsotropicElastic mat(200e3, 0.3, 3);
    
    EXPECT_DOUBLE_EQ(mat.getParameter("E"), 200e3);
    EXPECT_DOUBLE_EQ(mat.getParameter("nu"), 0.3);
    EXPECT_EQ(mat.strainSize(), 6);
}

TEST(IsotropicElasticTest, LameConstants) {
    IsotropicElastic mat(200e3, 0.3, 3);
    
    // λ = E*nu / ((1+nu)*(1-2*nu))
    Real lambda_expected = 200e3 * 0.3 / ((1.0 + 0.3) * (1.0 - 2.0 * 0.3));
    EXPECT_NEAR(mat.lambda(), lambda_expected, 1e-6);
    
    // μ = E / (2*(1+nu))
    Real mu_expected = 200e3 / (2.0 * (1.0 + 0.3));
    EXPECT_NEAR(mat.mu(), mu_expected, 1e-6);
}

TEST(IsotropicElasticTest, ElasticityTensor3D) {
    IsotropicElastic mat(200e3, 0.3, 3);
    DenseMatrix D = mat.elasticityTensor();
    
    EXPECT_EQ(D.rows(), 6);
    EXPECT_EQ(D.cols(), 6);
    
    // 检查对称性
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_NEAR(D(i, j), D(j, i), 1e-10);
        }
    }
    
    // 检查剪切模量（对角线后三个元素）
    Real mu = mat.mu();
    EXPECT_NEAR(D(3, 3), mu, 1e-6);
    EXPECT_NEAR(D(4, 4), mu, 1e-6);
    EXPECT_NEAR(D(5, 5), mu, 1e-6);
}

TEST(IsotropicElasticTest, UniaxialTension3D) {
    IsotropicElastic mat(200e3, 0.3, 3);
    StateVariables state = mat.createState();
    
    // 单轴拉伸：ε11 = 0.001
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.001;
    
    Vector stress(6, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 应力：σ11 ≠ 0, σ22 = σ33 ≈ 0 (小但非零，泊松效应)
    EXPECT_GT(stress[0], 0.0);
    
    // 检查应变能
    Vector total_strain = strain_inc;
    Real energy = mat.strainEnergy(total_strain, state);
    EXPECT_GT(energy, 0.0);
}

TEST(IsotropicElasticTest, PureShear2D) {
    IsotropicElastic mat(200e3, 0.3, 2, true);  // 平面应力
    StateVariables state = mat.createState();
    
    // 纯剪切：γ12 = 0.002
    Vector strain_inc(3, 0.0);
    strain_inc[2] = 0.002;
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 剪应力：σ12 = G * γ12
    Real G = mat.mu();
    EXPECT_NEAR(stress[2], G * 0.002, 1e-6);
    
    // 正应力应为零（纯剪切）
    EXPECT_NEAR(stress[0], 0.0, 1e-6);
    EXPECT_NEAR(stress[1], 0.0, 1e-6);
}

TEST(IsotropicElasticTest, InvalidParameters) {
    // 负杨氏模量
    EXPECT_THROW(IsotropicElastic(-200e3, 0.3, 3), std::invalid_argument);
    
    // 泊松比超出范围
    EXPECT_THROW(IsotropicElastic(200e3, 0.6, 3), std::invalid_argument);
    EXPECT_THROW(IsotropicElastic(200e3, -1.5, 3), std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════
// J2Plasticity 测试
// ═══════════════════════════════════════════════════════════════

TEST(J2PlasticityTest, Construction) {
    J2Plasticity mat(200e3, 0.3, 250.0, 1000.0, 3);
    
    EXPECT_DOUBLE_EQ(mat.getParameter("E"), 200e3);
    EXPECT_DOUBLE_EQ(mat.getParameter("nu"), 0.3);
    EXPECT_DOUBLE_EQ(mat.getParameter("sigma_y0"), 250.0);
    EXPECT_DOUBLE_EQ(mat.getParameter("H"), 1000.0);
}

TEST(J2PlasticityTest, YieldStress) {
    J2Plasticity mat(200e3, 0.3, 250.0, 1000.0, 3);
    
    // 初始屈服应力
    EXPECT_DOUBLE_EQ(mat.yieldStress(0.0), 250.0);
    
    // 硬化后
    Real eps_p = 0.01;
    Real sigma_y_expected = 250.0 + 1000.0 * eps_p;
    EXPECT_DOUBLE_EQ(mat.yieldStress(eps_p), sigma_y_expected);
}

TEST(J2PlasticityTest, VonMisesStress3D) {
    J2Plasticity mat(200e3, 0.3, 250.0, 0.0, 3);
    
    // 单轴应力状态：σ11 = 300
    Vector stress(6, 0.0);
    stress[0] = 300.0;
    
    Real q = mat.vonMisesStress(stress);
    
    // 单轴：q = |σ11|
    EXPECT_NEAR(q, 300.0, 1e-6);
}

TEST(J2PlasticityTest, ElasticLoading) {
    J2Plasticity mat(200e3, 0.3, 250.0, 0.0, 3);
    StateVariables state = mat.createState();
    
    // 小应变（弹性范围内）
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.0005;  // ε11 = 0.05% (小应变)
    
    Vector stress(6, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 应力应在弹性范围内
    Real q = mat.vonMisesStress(stress);
    EXPECT_LT(q, 250.0);  // 小于屈服应力
    
    // 塑性应变应为零
    EXPECT_DOUBLE_EQ(state.equiv_plastic_strain, 0.0);
}

TEST(J2PlasticityTest, PlasticLoading) {
    J2Plasticity mat(200e3, 0.3, 250.0, 1000.0, 3);
    StateVariables state = mat.createState();
    
    // 大应变（超过屈服）
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.002;  // ε11 = 0.2% (大应变)
    
    Vector stress(6, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 应力应在屈服面上
    Real q = mat.vonMisesStress(stress);
    Real sigma_y = mat.yieldStress(state.equiv_plastic_strain);
    
    // 屈服函数应接近零
    EXPECT_NEAR(q, sigma_y, 1.0);  // 允许小误差
    
    // 塑性应变应 > 0
    EXPECT_GT(state.equiv_plastic_strain, 0.0);
}

TEST(J2PlasticityTest, UnloadingBehavior) {
    J2Plasticity mat(200e3, 0.3, 250.0, 0.0, 3);
    StateVariables state = mat.createState();
    
    // 1. 加载到塑性
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.003;
    
    Vector stress(6, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    Real plastic_strain_after_loading = state.equiv_plastic_strain;
    EXPECT_GT(plastic_strain_after_loading, 0.0);
    
    // 2. 卸载（负应变增量）
    strain_inc[0] = -0.001;
    mat.computeStress(strain_inc, stress, state);
    
    // 塑性应变应保持不变（弹性卸载）
    EXPECT_DOUBLE_EQ(state.equiv_plastic_strain, plastic_strain_after_loading);
}

TEST(J2PlasticityTest, CyclicLoading) {
    J2Plasticity mat(200e3, 0.3, 250.0, 1000.0, 3);
    StateVariables state = mat.createState();
    
    Vector stress(6, 0.0);
    
    // 循环1：加载
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.002;
    mat.computeStress(strain_inc, stress, state);
    Real eps_p_1 = state.equiv_plastic_strain;
    
    // 循环2：卸载
    strain_inc[0] = -0.001;
    mat.computeStress(strain_inc, stress, state);
    EXPECT_DOUBLE_EQ(state.equiv_plastic_strain, eps_p_1);  // 弹性卸载
    
    // 循环3：再加载
    strain_inc[0] = 0.002;
    mat.computeStress(strain_inc, stress, state);
    
    // 塑性应变应累积
    EXPECT_GT(state.equiv_plastic_strain, eps_p_1);
}

TEST(J2PlasticityTest, InvalidParameters) {
    // 负屈服应力
    EXPECT_THROW(J2Plasticity(200e3, 0.3, -250.0, 0.0, 3), std::invalid_argument);
    
    // 负硬化模量
    EXPECT_THROW(J2Plasticity(200e3, 0.3, 250.0, -1000.0, 3), std::invalid_argument);
}

// main() provided by gtest_main
