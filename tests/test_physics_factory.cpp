#include <gtest/gtest.h>
#include "physics/physics_factory.h"
#include "material/isotropic_elastic.h"

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════
// PhysicsFactory 测试
// ═══════════════════════════════════════════════════════════

TEST(PhysicsFactoryTest, CreateHeat) {
    auto* physics = PhysicsFactory::create("Heat", {
        {"conductivity", 1.5},
        {"source", 10.0}
    });
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, CreateHeatWithDefaults) {
    auto* physics = PhysicsFactory::create("Heat");
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, CreateElasticity) {
    // 创建材料
    IsotropicElastic material(200e9, 0.3, 3);
    
    // 通过工厂创建弹性力学模块
    auto* physics = PhysicsFactory::create("Elasticity", 
        {{"dimension", 3}},
        {{"material_ptr", &material}}
    );
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, CreateElasticity2D) {
    IsotropicElastic material(70e9, 0.33, 2, true);
    
    auto* physics = PhysicsFactory::create("Elasticity",
        {{"dimension", 2}},
        {{"material_ptr", &material}}
    );
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, CreateElasticityConvenience) {
    IsotropicElastic material(200e9, 0.3, 3);
    
    auto* physics = PhysicsFactory::createElasticity(&material, 3);
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, CreateHeatConvenience) {
    auto* physics = PhysicsFactory::createHeat(2.5, 5.0);
    
    ASSERT_NE(physics, nullptr);
    delete physics;
}

TEST(PhysicsFactoryTest, MissingMaterialPointer) {
    EXPECT_THROW(
        PhysicsFactory::create("Elasticity", {{"dimension", 3}}),
        std::invalid_argument
    );
}

TEST(PhysicsFactoryTest, UnknownType) {
    EXPECT_THROW(
        PhysicsFactory::create("NonExistentPhysics"),
        std::invalid_argument
    );
}

TEST(PhysicsFactoryTest, IsRegistered) {
    EXPECT_TRUE(PhysicsFactory::isRegistered("Elasticity"));
    EXPECT_TRUE(PhysicsFactory::isRegistered("Heat"));
    EXPECT_FALSE(PhysicsFactory::isRegistered("UnknownPhysics"));
}

TEST(PhysicsFactoryTest, GetRegisteredTypes) {
    auto types = PhysicsFactory::getRegisteredTypes();
    
    EXPECT_GE(types.size(), 2);  // 至少 Elasticity, Heat
    
    bool found_elasticity = false;
    bool found_heat = false;
    
    for (const auto& name : types) {
        if (name == "Elasticity") found_elasticity = true;
        if (name == "Heat") found_heat = true;
    }
    
    EXPECT_TRUE(found_elasticity);
    EXPECT_TRUE(found_heat);
}

TEST(PhysicsFactoryTest, RegisterCustomPhysics) {
    // 注册一个自定义物理模块（实际上创建 HeatConductionUnified）
    PhysicsFactory::registerCreator("TestPhysics",
        [](const PhysicsFactory::Parameters& num_p,
           const PhysicsFactory::PointerParams& ptr_p) -> PhysicsBase* 
        {
            return new HeatConductionUnified(1.0, 0.0);
        }
    );
    
    EXPECT_TRUE(PhysicsFactory::isRegistered("TestPhysics"));
    
    auto* physics = PhysicsFactory::create("TestPhysics");
    ASSERT_NE(physics, nullptr);
    delete physics;
}
