#include <gtest/gtest.h>
#include "material/state_variables.h"
#include "material/material.h"
#include <sstream>

using namespace fem;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════════
// StateVariables 测试
// ═══════════════════════════════════════════════════════════════

TEST(StateVariablesTest, DefaultConstruction) {
    StateVariables sv;
    EXPECT_EQ(sv.plastic_strain.size(), 0);
    EXPECT_DOUBLE_EQ(sv.equiv_plastic_strain, 0.0);
    EXPECT_DOUBLE_EQ(sv.damage, 0.0);
}

TEST(StateVariablesTest, SizeConstruction) {
    StateVariables sv(6);  // 3D应力状态
    EXPECT_EQ(sv.plastic_strain.size(), 6);
    EXPECT_EQ(sv.back_stress.size(), 6);
    
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(sv.plastic_strain[i], 0.0);
        EXPECT_DOUBLE_EQ(sv.back_stress[i], 0.0);
    }
}

TEST(StateVariablesTest, ScalarVariables) {
    StateVariables sv;
    
    sv.setScalar("kappa", 100.0);
    sv.setScalar("temperature", 300.0);
    
    EXPECT_DOUBLE_EQ(sv.getScalar("kappa"), 100.0);
    EXPECT_DOUBLE_EQ(sv.getScalar("temperature"), 300.0);
    EXPECT_DOUBLE_EQ(sv.getScalar("nonexistent", 42.0), 42.0);  // 默认值
}

TEST(StateVariablesTest, TensorVariables) {
    StateVariables sv;
    
    Vector alpha(6, 1.5);
    sv.setTensor("kinematic_hardening", alpha);
    
    EXPECT_TRUE(sv.hasTensor("kinematic_hardening"));
    EXPECT_FALSE(sv.hasTensor("nonexistent"));
    
    Vector retrieved = sv.getTensor("kinematic_hardening");
    EXPECT_EQ(retrieved.size(), 6);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(retrieved[i], 1.5);
    }
}

TEST(StateVariablesTest, Reset) {
    StateVariables sv(6);
    sv.equiv_plastic_strain = 0.5;
    sv.damage = 0.2;
    sv.plastic_strain[0] = 0.01;
    sv.setScalar("test", 123.0);
    
    sv.reset();
    
    EXPECT_DOUBLE_EQ(sv.equiv_plastic_strain, 0.0);
    EXPECT_DOUBLE_EQ(sv.damage, 0.0);
    EXPECT_DOUBLE_EQ(sv.plastic_strain[0], 0.0);
    EXPECT_TRUE(sv.scalar_vars.empty());
}

TEST(StateVariablesTest, Clone) {
    StateVariables sv(6);
    sv.equiv_plastic_strain = 0.1;
    sv.damage = 0.05;
    sv.plastic_strain[2] = 0.002;
    sv.setScalar("hardening", 50.0);
    
    StateVariables copy = sv.clone();
    
    EXPECT_DOUBLE_EQ(copy.equiv_plastic_strain, 0.1);
    EXPECT_DOUBLE_EQ(copy.damage, 0.05);
    EXPECT_DOUBLE_EQ(copy.plastic_strain[2], 0.002);
    EXPECT_DOUBLE_EQ(copy.getScalar("hardening"), 50.0);
    
    // 验证深拷贝：修改原始不影响副本
    sv.damage = 0.9;
    EXPECT_DOUBLE_EQ(copy.damage, 0.05);
}

TEST(StateVariablesTest, Serialization) {
    StateVariables sv(6);
    sv.equiv_plastic_strain = 0.123;
    sv.damage = 0.456;
    sv.plastic_strain[0] = 0.001;
    sv.plastic_strain[3] = 0.002;
    sv.setScalar("kappa", 99.9);
    
    Vector tensor(3, 7.7);
    sv.setTensor("alpha", tensor);
    
    // 序列化
    std::ostringstream oss(std::ios::binary);
    sv.serialize(oss);
    
    // 反序列化
    StateVariables loaded;
    std::istringstream iss(oss.str(), std::ios::binary);
    loaded.deserialize(iss);
    
    // 验证
    EXPECT_DOUBLE_EQ(loaded.equiv_plastic_strain, 0.123);
    EXPECT_DOUBLE_EQ(loaded.damage, 0.456);
    EXPECT_EQ(loaded.plastic_strain.size(), 6);
    EXPECT_DOUBLE_EQ(loaded.plastic_strain[0], 0.001);
    EXPECT_DOUBLE_EQ(loaded.plastic_strain[3], 0.002);
    EXPECT_DOUBLE_EQ(loaded.getScalar("kappa"), 99.9);
    
    Vector loaded_tensor = loaded.getTensor("alpha");
    EXPECT_EQ(loaded_tensor.size(), 3);
    EXPECT_DOUBLE_EQ(loaded_tensor[0], 7.7);
}

// ═══════════════════════════════════════════════════════════════
// Material 基类测试 (使用 Mock 子类)
// ═══════════════════════════════════════════════════════════════

class MockMaterial : public Material {
public:
    MockMaterial() : Material(6) {
        setParameter("E", 200e3);
        setParameter("nu", 0.3);
    }
    
    void computeStress(const Vector&, Vector&, StateVariables&) override {
        // Mock: 不实现
    }
    
    void computeTangent(const Vector&, DenseMatrix&, const StateVariables&) override {
        // Mock: 不实现
    }
    
    Real strainEnergy(const Vector&, const StateVariables&) const override {
        return 0.0;
    }
    
    StateVariables createState() const override {
        return StateVariables(6);
    }
    
    std::string typeName() const override {
        return "MockMaterial";
    }
};

TEST(MaterialTest, ParameterManagement) {
    MockMaterial mat;
    
    EXPECT_DOUBLE_EQ(mat.getParameter("E"), 200e3);
    EXPECT_DOUBLE_EQ(mat.getParameter("nu"), 0.3);
    
    mat.setParameter("sigma_y", 250.0);
    EXPECT_DOUBLE_EQ(mat.getParameter("sigma_y"), 250.0);
    
    auto names = mat.getParameterNames();
    EXPECT_EQ(names.size(), 3);  // E, nu, sigma_y
    EXPECT_EQ(names[0], "E");
    EXPECT_EQ(names[1], "nu");
    EXPECT_EQ(names[2], "sigma_y");
}

TEST(MaterialTest, ParameterNotFound) {
    MockMaterial mat;
    EXPECT_THROW(mat.getParameter("nonexistent"), std::invalid_argument);
}

TEST(MaterialTest, StrainSize) {
    MockMaterial mat;
    EXPECT_EQ(mat.strainSize(), 6);  // 3D
}

TEST(MaterialTest, TypeName) {
    MockMaterial mat;
    EXPECT_EQ(mat.typeName(), "MockMaterial");
}

TEST(MaterialTest, CreateState) {
    MockMaterial mat;
    StateVariables sv = mat.createState();
    EXPECT_EQ(sv.plastic_strain.size(), 6);
    EXPECT_DOUBLE_EQ(sv.equiv_plastic_strain, 0.0);
}

TEST(MaterialTest, GeometricStiffnessDefault) {
    MockMaterial mat;
    Vector stress(6, 100.0);
    DenseMatrix K_geo;
    
    mat.computeGeometricStiffness(stress, K_geo);
    
    EXPECT_EQ(K_geo.rows(), 6);
    EXPECT_EQ(K_geo.cols(), 6);
    
    // 默认实现：零矩阵
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_DOUBLE_EQ(K_geo(i, j), 0.0);
        }
    }
}

// main() provided by gtest_main
