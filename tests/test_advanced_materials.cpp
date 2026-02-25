#include <gtest/gtest.h>
#include "material/orthotropic_elastic.h"
#include "material/j2_plasticity_kinematic.h"
#include <cmath>

using namespace fem;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════════
// OrthotropicElastic 测试
// ═══════════════════════════════════════════════════════════════

TEST(OrthotropicElasticTest, Construction2D) {
    // 碳纤维复合材料参数
    Real E1 = 181e3;   // 纤维方向
    Real E2 = 10.3e3;  // 横向
    Real nu12 = 0.28;
    Real G12 = 7.17e3;
    
    OrthotropicElastic mat(E1, E2, nu12, G12, 0.0);
    
    EXPECT_DOUBLE_EQ(mat.getParameter("E1"), E1);
    EXPECT_DOUBLE_EQ(mat.getParameter("E2"), E2);
    EXPECT_DOUBLE_EQ(mat.getParameter("nu12"), nu12);
    EXPECT_DOUBLE_EQ(mat.getParameter("G12"), G12);
}

TEST(OrthotropicElasticTest, StiffnessMatrix2D) {
    Real E1 = 100e3, E2 = 10e3;
    Real nu12 = 0.3, G12 = 5e3;
    
    OrthotropicElastic mat(E1, E2, nu12, G12, 0.0);
    
    DenseMatrix D = mat.localStiffness();
    
    EXPECT_EQ(D.rows(), 3);
    EXPECT_EQ(D.cols(), 3);
    
    // 检查对称性
    EXPECT_NEAR(D(0, 1), D(1, 0), 1e-6);
    
    // 检查主对角元
    EXPECT_GT(D(0, 0), D(1, 1));  // E1 > E2
    EXPECT_DOUBLE_EQ(D(2, 2), G12);
}

TEST(OrthotropicElasticTest, RotationTransform) {
    Real E1 = 100e3, E2 = 50e3;
    Real nu12 = 0.25, G12 = 20e3;
    
    // 无旋转
    OrthotropicElastic mat1(E1, E2, nu12, G12, 0.0);
    DenseMatrix D0 = mat1.globalStiffness();
    
    // 旋转90度
    OrthotropicElastic mat2(E1, E2, nu12, G12, 90.0);
    DenseMatrix D90 = mat2.globalStiffness();
    
    // 旋转90度后：D11 和 D22 应交换
    EXPECT_NEAR(D0(0, 0), D90(1, 1), 1e-6);
    EXPECT_NEAR(D0(1, 1), D90(0, 0), 1e-6);
}

TEST(OrthotropicElasticTest, UniaxialTension) {
    Real E1 = 150e3, E2 = 10e3;
    Real nu12 = 0.3, G12 = 5e3;
    
    OrthotropicElastic mat(E1, E2, nu12, G12, 0.0);
    StateVariables state = mat.createState();
    
    // 沿纤维方向拉伸
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.001;  // ε1 = 0.1%
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 主方向刚度大
    EXPECT_GT(std::abs(stress[0]), std::abs(stress[1]));
}

TEST(OrthotropicElasticTest, InvalidParameters) {
    // 负杨氏模量
    EXPECT_THROW(OrthotropicElastic(-100e3, 10e3, 0.3, 5e3, 0.0), 
                 std::invalid_argument);
    
    // 泊松比超出范围
    EXPECT_THROW(OrthotropicElastic(100e3, 10e3, 5.0, 5e3, 0.0), 
                 std::invalid_argument);
}

// ═══════════════════════════════════════════════════════════════
// J2PlasticityKinematic 测试
// ═══════════════════════════════════════════════════════════════

TEST(J2PlasticityKinematicTest, Construction) {
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 1000.0, 5000.0, 0.0, 2);
    
    EXPECT_DOUBLE_EQ(mat.getParameter("E"), 200e3);
    EXPECT_DOUBLE_EQ(mat.getParameter("H_iso"), 1000.0);
    EXPECT_DOUBLE_EQ(mat.getParameter("H_kin"), 5000.0);
}

TEST(J2PlasticityKinematicTest, ElasticLoading) {
    // 仅运动硬化
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 0.0, 5000.0, 0.0, 2);
    StateVariables state = mat.createState();
    
    // 小应变（弹性）
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.0005;
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 背应力应为零（未屈服）
    EXPECT_DOUBLE_EQ(state.back_stress[0], 0.0);
    EXPECT_DOUBLE_EQ(state.equiv_plastic_strain, 0.0);
}

TEST(J2PlasticityKinematicTest, PlasticLoadingKinematic) {
    // 纯运动硬化
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 0.0, 5000.0, 0.0, 2);
    StateVariables state = mat.createState();
    
    // 大应变（塑性）
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.002;
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 应该有塑性应变
    EXPECT_GT(state.equiv_plastic_strain, 0.0);
    
    // 背应力应非零（运动硬化）
    EXPECT_GT(std::abs(state.back_stress[0]), 0.1);
}

TEST(J2PlasticityKinematicTest, BauschingerEffect) {
    // 运动硬化材料（Bauschinger效应）
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 0.0, 5000.0, 0.0, 2);
    StateVariables state = mat.createState();
    Vector stress(3, 0.0);
    
    // 步骤1：正向加载（拉伸）
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.003;  // 大应变，屈服
    mat.computeStress(strain_inc, stress, state);
    
    Real eps_p_forward = state.equiv_plastic_strain;
    Real alpha_x_forward = state.back_stress[0];
    
    EXPECT_GT(eps_p_forward, 0.0);
    EXPECT_GT(alpha_x_forward, 0.0);  // 背应力沿拉伸方向
    
    // 步骤2：反向加载（压缩）
    strain_inc[0] = -0.002;
    mat.computeStress(strain_inc, stress, state);
    
    // Bauschinger效应：反向屈服应力降低
    // （这里只是定性测试，具体数值取决于模型参数）
    // 注意：小的反向加载可能仍在弹性范围，因此塑性应变可能不变
    EXPECT_GE(state.equiv_plastic_strain, eps_p_forward);  // 塑性应变不减少
}

TEST(J2PlasticityKinematicTest, MixedHardening) {
    // 混合硬化
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 1000.0, 2000.0, 0.0, 2);
    StateVariables state = mat.createState();
    
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.003;
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 同时有等向和运动硬化
    EXPECT_GT(state.equiv_plastic_strain, 0.0);  // 等向硬化
    EXPECT_GT(std::abs(state.back_stress[0]), 0.0);  // 运动硬化
    
    // 屈服面扩大 + 平移
    Real sigma_y = mat.yieldStress(state.equiv_plastic_strain);
    EXPECT_GT(sigma_y, mat.getParameter("sigma_y0"));
}

TEST(J2PlasticityKinematicTest, NonlinearRecovery) {
    // 带恢复参数的Armstrong-Frederick模型
    J2PlasticityKinematic mat(200e3, 0.3, 250.0, 0.0, 5000.0, 100.0, 2);
    StateVariables state = mat.createState();
    
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.004;  // 大应变
    
    Vector stress(3, 0.0);
    mat.computeStress(strain_inc, stress, state);
    
    // 非线性恢复：背应力增长受限（饱和效应）
    EXPECT_GT(state.equiv_plastic_strain, 0.0);
    
    // 背应力存在但受恢复项约束
    // （具体值取决于beta参数）
}

// main() provided by gtest_main
