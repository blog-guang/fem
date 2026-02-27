/**
 * test_data.cpp - 数据管理系统单元测试
 */

#include <gtest/gtest.h>
#include "data/field_data.h"
#include "data/data_manager.h"

using namespace fem;
using namespace fem::data;

// ═══════════════════════════════════════════════════════════════
// RealData Tests
// ═══════════════════════════════════════════════════════════════

TEST(RealDataTest, Construction) {
    RealData field("temperature", DataLocation::Node, 10, 20.0);
    
    EXPECT_EQ(field.name(), "temperature");
    EXPECT_EQ(field.location(), DataLocation::Node);
    EXPECT_EQ(field.size(), 10);
    EXPECT_EQ(field.type_name(), "Real");
    EXPECT_EQ(field.default_value(), 20.0);
}

TEST(RealDataTest, SetGet) {
    RealData field("temp", DataLocation::Node, 100);
    
    field.set(10, 100.5);
    EXPECT_DOUBLE_EQ(field.get(10), 100.5);
    EXPECT_DOUBLE_EQ(field[10], 100.5);
}

TEST(RealDataTest, Resize) {
    RealData field("temp", DataLocation::Node, 10);
    EXPECT_EQ(field.size(), 10);
    
    field.resize(50);
    EXPECT_EQ(field.size(), 50);
}

TEST(RealDataTest, Reset) {
    RealData field("temp", DataLocation::Node, 10, 20.0);
    
    field.set(5, 100.0);
    EXPECT_DOUBLE_EQ(field[5], 100.0);
    
    field.reset();
    EXPECT_DOUBLE_EQ(field[5], 20.0);
}

TEST(RealDataTest, GaussPoint) {
    Index n_elems = 10;
    Index n_gp = 4;
    
    RealData field("damage", DataLocation::GaussPoint, n_elems * n_gp);
    
    field.set_gauss_point(5, 2, n_gp, 0.8);
    EXPECT_DOUBLE_EQ(field.get_gauss_point(5, 2, n_gp), 0.8);
}

// ═══════════════════════════════════════════════════════════════
// VectorData Tests
// ═══════════════════════════════════════════════════════════════

TEST(VectorDataTest, Construction) {
    VectorData field("displacement", DataLocation::Node, 10);
    
    EXPECT_EQ(field.name(), "displacement");
    EXPECT_EQ(field.type_name(), "Vector");
    EXPECT_EQ(field.size(), 10);
}

TEST(VectorDataTest, SetGet) {
    VectorData field("disp", DataLocation::Node, 10);
    
    Vector u(3);
    u[0] = 1.0;
    u[1] = 2.0;
    u[2] = 3.0;
    
    field.set(5, u);
    
    const Vector& u_read = field.get(5);
    EXPECT_DOUBLE_EQ(u_read[0], 1.0);
    EXPECT_DOUBLE_EQ(u_read[1], 2.0);
    EXPECT_DOUBLE_EQ(u_read[2], 3.0);
}

// ═══════════════════════════════════════════════════════════════
// IntData Tests
// ═══════════════════════════════════════════════════════════════

TEST(IntDataTest, Construction) {
    IntData field("material_id", DataLocation::Element, 20, 1);
    
    EXPECT_EQ(field.type_name(), "Int");
    EXPECT_EQ(field.default_value(), 1);
}

TEST(IntDataTest, SetGet) {
    IntData field("mat_id", DataLocation::Element, 10);
    
    field.set(3, 5);
    EXPECT_EQ(field.get(3), 5);
}

// ═══════════════════════════════════════════════════════════════
// BoolData Tests
// ═══════════════════════════════════════════════════════════════

TEST(BoolDataTest, Construction) {
    BoolData field("is_active", DataLocation::Element, 10, true);
    
    EXPECT_EQ(field.type_name(), "Bool");
    EXPECT_TRUE(field.default_value());
}

TEST(BoolDataTest, SetGet) {
    BoolData field("active", DataLocation::Element, 10, true);
    
    field.set(5, false);
    EXPECT_FALSE(field.get(5));
    EXPECT_FALSE(field[5]);
    EXPECT_TRUE(field.get(3));
    EXPECT_TRUE(field[3]);
}

// ═══════════════════════════════════════════════════════════════
// DataManager Tests
// ═══════════════════════════════════════════════════════════════

TEST(DataManagerTest, AddField) {
    DataManager manager;
    
    auto* temp = manager.add_field<RealData>("temperature", DataLocation::Node, 100);
    
    ASSERT_NE(temp, nullptr);
    EXPECT_EQ(temp->name(), "temperature");
    EXPECT_TRUE(manager.has_field("temperature"));
    EXPECT_EQ(manager.num_fields(), 1);
}

TEST(DataManagerTest, GetField) {
    DataManager manager;
    
    manager.add_field<RealData>("temp", DataLocation::Node, 10);
    
    auto* temp = manager.get_field<RealData>("temp");
    ASSERT_NE(temp, nullptr);
    
    temp->set(5, 100.0);
    EXPECT_DOUBLE_EQ(temp->get(5), 100.0);
}

TEST(DataManagerTest, RemoveField) {
    DataManager manager;
    
    manager.add_field<RealData>("temp", DataLocation::Node, 10);
    EXPECT_TRUE(manager.has_field("temp"));
    
    manager.remove_field("temp");
    EXPECT_FALSE(manager.has_field("temp"));
    EXPECT_EQ(manager.num_fields(), 0);
}

TEST(DataManagerTest, MultipleFields) {
    DataManager manager;
    
    manager.add_field<RealData>("temperature", DataLocation::Node, 100);
    manager.add_field<VectorData>("displacement", DataLocation::Node, 100);
    manager.add_field<IntData>("material_id", DataLocation::Element, 80);
    
    EXPECT_EQ(manager.num_fields(), 3);
    EXPECT_EQ(manager.num_fields(DataLocation::Node), 2);
    EXPECT_EQ(manager.num_fields(DataLocation::Element), 1);
}

TEST(DataManagerTest, ResizeAll) {
    DataManager manager;
    
    manager.add_field<RealData>("temp", DataLocation::Node, 100);
    manager.add_field<VectorData>("disp", DataLocation::Node, 100);
    manager.add_field<IntData>("mat_id", DataLocation::Element, 80);
    
    manager.resize_all(DataLocation::Node, 200);
    
    auto* temp = manager.get_field<RealData>("temp");
    auto* disp = manager.get_field<VectorData>("disp");
    auto* mat_id = manager.get_field<IntData>("mat_id");
    
    EXPECT_EQ(temp->size(), 200);
    EXPECT_EQ(disp->size(), 200);
    EXPECT_EQ(mat_id->size(), 80);  // 不变
}

TEST(DataManagerTest, GetFieldNames) {
    DataManager manager;
    
    manager.add_field<RealData>("temp", DataLocation::Node, 10);
    manager.add_field<VectorData>("disp", DataLocation::Node, 10);
    manager.add_field<IntData>("mat", DataLocation::Element, 5);
    
    auto names = manager.get_field_names();
    EXPECT_EQ(names.size(), 3);
    
    auto node_names = manager.get_field_names(DataLocation::Node);
    EXPECT_EQ(node_names.size(), 2);
}

TEST(DataManagerTest, ClearAll) {
    DataManager manager;
    
    manager.add_field<RealData>("temp", DataLocation::Node, 100);
    manager.add_field<VectorData>("disp", DataLocation::Node, 100);
    
    manager.clear_all();
    
    auto* temp = manager.get_field<RealData>("temp");
    EXPECT_EQ(temp->size(), 0);
}

TEST(DataManagerTest, ResetAll) {
    DataManager manager;
    
    auto* temp = manager.add_field<RealData>("temp", DataLocation::Node, 10, 20.0);
    temp->set(5, 100.0);
    
    manager.reset_all();
    EXPECT_DOUBLE_EQ(temp->get(5), 20.0);
}
