#include <gtest/gtest.h>
#include "physics/elasticity_nonlinear.h"
#include "material/isotropic_elastic.h"
#include "material/neo_hookean.h"

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════
// ElasticityNonlinear 基础测试
// ═══════════════════════════════════════════════════════════

TEST(ElasticityNonlinearTest, Construction) {
    // 创建材料
    IsotropicElastic mat(200e9, 0.3, 3, false);
    
    // 创建几何非线性物理模块
    ElasticityNonlinear physics(&mat, 3, true);
    
    EXPECT_EQ(physics.dimension(), 3);
    EXPECT_TRUE(physics.is_geometric_stiffness_enabled());
}

TEST(ElasticityNonlinearTest, Construction_2D) {
    IsotropicElastic mat(200e9, 0.3, 2, false);
    ElasticityNonlinear physics(&mat, 2, true);
    
    EXPECT_EQ(physics.dimension(), 2);
}

TEST(ElasticityNonlinearTest, DisableGeometricStiffness) {
    IsotropicElastic mat(200e9, 0.3, 3, false);
    ElasticityNonlinear physics(&mat, 3, false);
    
    EXPECT_FALSE(physics.is_geometric_stiffness_enabled());
    
    physics.set_geometric_stiffness(true);
    EXPECT_TRUE(physics.is_geometric_stiffness_enabled());
}

TEST(ElasticityNonlinearTest, StressUpdateMethod) {
    IsotropicElastic mat(200e9, 0.3, 3, false);
    ElasticityNonlinear physics(&mat, 3, false);
    
    // 设置应力更新方法
    physics.set_stress_update_method(
        ElasticityNonlinear::StressUpdate::GREEN_LAGRANGE
    );
    
    // 验证设置成功（通过能够调用即可）
    SUCCEED();
}

TEST(ElasticityNonlinearTest, InvalidDimension) {
    IsotropicElastic mat(200e9, 0.3, 3, false);
    
    // 维度必须是 2 或 3
    EXPECT_THROW(
        ElasticityNonlinear(&mat, 1, false),
        std::invalid_argument
    );
    
    EXPECT_THROW(
        ElasticityNonlinear(&mat, 4, false),
        std::invalid_argument
    );
}

TEST(ElasticityNonlinearTest, NullMaterial) {
    // 材料指针不能为空
    EXPECT_THROW(
        ElasticityNonlinear(nullptr, 3, false),
        std::invalid_argument
    );
}

TEST(ElasticityNonlinearTest, NeoHookeanMaterial) {
    // 使用 Neo-Hookean 材料
    NeoHookean mat(10e6, 0.3, 3, true);
    ElasticityNonlinear physics(&mat, 3, false);
    
    EXPECT_EQ(physics.dimension(), 3);
    EXPECT_NE(physics.material(), nullptr);
}

TEST(ElasticityNonlinearTest, MaterialPointer) {
    IsotropicElastic mat(200e9, 0.3, 3, false);
    ElasticityNonlinear physics(&mat, 3, false);
    
    // 验证材料指针
    EXPECT_EQ(physics.material(), &mat);
}
