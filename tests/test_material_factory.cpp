#include <gtest/gtest.h>
#include "material/material_factory.h"
#include "material/isotropic_elastic.h"
#include "material/j2_plasticity.h"

using namespace fem;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════
// MaterialFactory 测试
// ═══════════════════════════════════════════════════════════

TEST(MaterialFactoryTest, CreateIsotropicElastic) {
    auto mat = MaterialFactory::create("IsotropicElastic", {
        {"E", 200e9},
        {"nu", 0.3},
        {"dimension", 3}
    });
    
    ASSERT_NE(mat, nullptr);
    // typeName() may include dimension info
    EXPECT_TRUE(mat->typeName().find("IsotropicElastic") != std::string::npos);
    EXPECT_EQ(mat->strainSize(), 6);  // 3D
}

TEST(MaterialFactoryTest, CreateIsotropicElastic2D) {
    auto mat = MaterialFactory::create("IsotropicElastic", {
        {"E", 70e9},
        {"nu", 0.33},
        {"dimension", 2},
        {"plane_stress", 1}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_EQ(mat->strainSize(), 3);  // 2D
}

TEST(MaterialFactoryTest, CreateJ2Plasticity) {
    auto mat = MaterialFactory::create("J2Plasticity", {
        {"E", 200e9},
        {"nu", 0.3},
        {"sigma_y", 250e6},
        {"H", 1e9},
        {"dimension", 3}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_TRUE(mat->typeName().find("J2Plasticity") != std::string::npos);
}

TEST(MaterialFactoryTest, MissingParameter) {
    EXPECT_THROW(
        MaterialFactory::create("IsotropicElastic", {
            {"E", 200e9}
            // 缺少 nu
        }),
        std::invalid_argument
    );
}

TEST(MaterialFactoryTest, UnknownType) {
    EXPECT_THROW(
        MaterialFactory::create("NonExistentMaterial", {
            {"param", 1.0}
        }),
        std::invalid_argument
    );
}

TEST(MaterialFactoryTest, IsRegistered) {
    EXPECT_TRUE(MaterialFactory::isRegistered("IsotropicElastic"));
    EXPECT_TRUE(MaterialFactory::isRegistered("J2Plasticity"));
    EXPECT_FALSE(MaterialFactory::isRegistered("UnknownMaterial"));
}

TEST(MaterialFactoryTest, GetRegisteredTypes) {
    auto types = MaterialFactory::getRegisteredTypes();
    
    EXPECT_GE(types.size(), 3);  // 至少有 IsotropicElastic, J2Plasticity 等
    
    bool found_elastic = false;
    bool found_plastic = false;
    
    for (const auto& name : types) {
        if (name == "IsotropicElastic") found_elastic = true;
        if (name == "J2Plasticity") found_plastic = true;
    }
    
    EXPECT_TRUE(found_elastic);
    EXPECT_TRUE(found_plastic);
}

TEST(MaterialFactoryTest, DefaultParameters) {
    // H 默认为 0（理想塑性）
    auto mat = MaterialFactory::create("J2Plasticity", {
        {"E", 200e9},
        {"nu", 0.3},
        {"sigma_y", 250e6}
        // H 使用默认值 0
    });
    
    ASSERT_NE(mat, nullptr);
}

TEST(MaterialFactoryTest, RegisterCustomMaterial) {
    // 注册一个简单的自定义材料创建器
    MaterialFactory::registerCreator("TestMaterial", 
        [](const MaterialFactory::Parameters& p) -> MaterialPtr {
            // 创建一个 IsotropicElastic 作为测试
            Real E = p.at("E");
            Real nu = p.at("nu");
            return std::make_shared<IsotropicElastic>(E, nu, 3);
        }
    );
    
    EXPECT_TRUE(MaterialFactory::isRegistered("TestMaterial"));
    
    auto mat = MaterialFactory::create("TestMaterial", {
        {"E", 100e9},
        {"nu", 0.25}
    });
    
    ASSERT_NE(mat, nullptr);
}

TEST(MaterialFactoryTest, MaterialStressComputation) {
    // 通过工厂创建材料，测试功能正常
    auto mat = MaterialFactory::create("IsotropicElastic", {
        {"E", 200e9},
        {"nu", 0.3},
        {"dimension", 3}
    });
    
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.001;  // ε11 拉伸
    
    Vector stress(6);
    StateVariables state = mat->createState();
    
    mat->computeStress(strain_inc, stress, state);
    
    EXPECT_GT(stress[0], 0.0);  // σ11 > 0（主向拉伸应力）
    EXPECT_GT(stress[1], 0.0);  // σ22 > 0（3D情况下，由于泊松约束）
    EXPECT_GT(stress[2], 0.0);  // σ33 > 0（同上）
}

TEST(MaterialFactoryTest, OrthotropicElastic2D) {
    auto mat = MaterialFactory::create("OrthotropicElastic", {
        {"E1", 150e9},
        {"E2", 100e9},
        {"nu12", 0.25},
        {"G12", 50e9},
        {"dimension", 2}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_TRUE(mat->typeName().find("OrthotropicElastic") != std::string::npos);
    EXPECT_EQ(mat->strainSize(), 3);
}

TEST(MaterialFactoryTest, J2PlasticityKinematic) {
    auto mat = MaterialFactory::create("J2PlasticityKinematic", {
        {"E", 200e9},
        {"nu", 0.3},
        {"sigma_y", 250e6},
        {"H_iso", 1e9},
        {"H_kin", 2e9},
        {"dimension", 3}
    });
    
    ASSERT_NE(mat, nullptr);
    EXPECT_TRUE(mat->typeName().find("J2PlasticityKinematic") != std::string::npos);
}
