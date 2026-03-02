#include <gtest/gtest.h>
#include "assembly/load_step.h"

using namespace fem;

// ═══════════════════════════════════════════════════════════
// LoadStep 测试
// ═══════════════════════════════════════════════════════════

TEST(LoadStepTest, Construction) {
    LoadStep step(1);
    
    EXPECT_EQ(step.id(), 1);
    EXPECT_EQ(step.num_substeps(), 1);
    EXPECT_EQ(step.ramp_type(), LoadRampType::RAMPED);
    EXPECT_DOUBLE_EQ(step.time_start(), 0.0);
    EXPECT_DOUBLE_EQ(step.time_end(), 1.0);
}

TEST(LoadStepTest, SetParameters) {
    LoadStep step;
    
    step.set_id(2);
    step.set_num_substeps(10);
    step.set_ramp_type(LoadRampType::STEPPED);
    step.set_time(0.0, 2.0);
    
    EXPECT_EQ(step.id(), 2);
    EXPECT_EQ(step.num_substeps(), 10);
    EXPECT_EQ(step.ramp_type(), LoadRampType::STEPPED);
    EXPECT_DOUBLE_EQ(step.time_start(), 0.0);
    EXPECT_DOUBLE_EQ(step.time_end(), 2.0);
}

TEST(LoadStepTest, DisplacementBC) {
    LoadStep step;
    
    // 添加位移边界条件：节点 5, DOF X, 0 → 0.01 m
    step.set_displacement(5, 0, 0.0, 0.01);
    
    // 在 t=0 时刻
    auto bcs_start = step.get_displacements(0.0);
    auto key = std::make_pair(Index(5), 0);
    EXPECT_DOUBLE_EQ(bcs_start[key], 0.0);
    
    // 在 t=0.5 时刻（中间）
    auto bcs_mid = step.get_displacements(0.5);
    EXPECT_DOUBLE_EQ(bcs_mid[key], 0.005);
    
    // 在 t=1.0 时刻
    auto bcs_end = step.get_displacements(1.0);
    EXPECT_DOUBLE_EQ(bcs_end[key], 0.01);
}

TEST(LoadStepTest, ForceLoad) {
    LoadStep step;
    
    // 添加力载荷：节点 10, DOF Z, 0 → 1000 N
    step.set_force(10, 2, 0.0, 1000.0);
    
    auto forces = step.get_forces(0.5);
    auto key = std::make_pair(Index(10), 2);
    EXPECT_DOUBLE_EQ(forces[key], 500.0);
}

TEST(LoadStepTest, RampedInterpolation) {
    LoadStep step;
    step.set_num_substeps(10);
    step.set_ramp_type(LoadRampType::RAMPED);
    step.set_displacement(1, 0, 0.0, 1.0);
    
    auto key = std::make_pair(Index(1), 0);
    
    // 子步 1: t=0.1
    auto bcs_1 = step.get_displacements(0.1);
    EXPECT_NEAR(bcs_1[key], 0.1, 1e-10);
    
    // 子步 5: t=0.5
    auto bcs_5 = step.get_displacements(0.5);
    EXPECT_NEAR(bcs_5[key], 0.5, 1e-10);
    
    // 子步 10: t=1.0
    auto bcs_10 = step.get_displacements(1.0);
    EXPECT_NEAR(bcs_10[key], 1.0, 1e-10);
}

TEST(LoadStepTest, SteppedInterpolation) {
    LoadStep step;
    step.set_ramp_type(LoadRampType::STEPPED);
    step.set_displacement(1, 0, 0.0, 1.0);
    
    auto key = std::make_pair(Index(1), 0);
    
    // 阶跃：t>0 时立即跳到终值
    auto bcs_start = step.get_displacements(0.0);
    EXPECT_DOUBLE_EQ(bcs_start[key], 0.0);
    
    auto bcs_epsilon = step.get_displacements(0.01);
    EXPECT_DOUBLE_EQ(bcs_epsilon[key], 1.0);
    
    auto bcs_mid = step.get_displacements(0.5);
    EXPECT_DOUBLE_EQ(bcs_mid[key], 1.0);
}

TEST(LoadStepTest, SubstepTime) {
    LoadStep step;
    step.set_num_substeps(4);
    step.set_time(0.0, 1.0);
    
    EXPECT_DOUBLE_EQ(step.substep_time(1), 0.25);
    EXPECT_DOUBLE_EQ(step.substep_time(2), 0.50);
    EXPECT_DOUBLE_EQ(step.substep_time(3), 0.75);
    EXPECT_DOUBLE_EQ(step.substep_time(4), 1.00);
}

TEST(LoadStepTest, LoadScale) {
    LoadStep step;
    step.set_displacement(1, 0, 0.0, 1.0);
    step.set_load_scale(2.0);  // 载荷缩放因子 = 2
    
    auto bcs = step.get_displacements(1.0);
    auto key = std::make_pair(Index(1), 0);
    EXPECT_DOUBLE_EQ(bcs[key], 2.0);  // 1.0 * 2.0
}

TEST(LoadStepTest, MultipleBC) {
    LoadStep step;
    
    step.set_displacement(1, 0, 0.0, 0.01);  // 节点 1, X
    step.set_displacement(1, 1, 0.0, 0.02);  // 节点 1, Y
    step.set_displacement(2, 2, 0.0, 0.03);  // 节点 2, Z
    
    auto bcs = step.get_displacements(0.5);
    
    auto key1 = std::make_pair(Index(1), 0);
    auto key2 = std::make_pair(Index(1), 1);
    auto key3 = std::make_pair(Index(2), 2);
    
    EXPECT_DOUBLE_EQ(bcs[key1], 0.005);
    EXPECT_DOUBLE_EQ(bcs[key2], 0.010);
    EXPECT_DOUBLE_EQ(bcs[key3], 0.015);
}

TEST(LoadStepTest, RemoveBC) {
    LoadStep step;
    
    step.set_displacement(1, 0, 0.0, 1.0);
    step.set_displacement(2, 1, 0.0, 2.0);
    
    step.remove_displacement(1, 0);
    
    auto bcs = step.get_displacements(1.0);
    
    auto key1 = std::make_pair(Index(1), 0);
    auto key2 = std::make_pair(Index(2), 1);
    
    EXPECT_EQ(bcs.count(key1), 0);  // 已删除
    EXPECT_EQ(bcs.count(key2), 1);  // 仍存在
}

// ═══════════════════════════════════════════════════════════
// LoadStepManager 测试
// ═══════════════════════════════════════════════════════════

TEST(LoadStepManagerTest, AddLoadStep) {
    LoadStepManager manager;
    
    LoadStep step1(1);
    LoadStep step2(2);
    
    manager.add_load_step(step1);
    manager.add_load_step(step2);
    
    EXPECT_EQ(manager.num_load_steps(), 2);
}

TEST(LoadStepManagerTest, GetLoadStep) {
    LoadStepManager manager;
    
    LoadStep step1(1);
    step1.set_num_substeps(5);
    
    LoadStep step2(2);
    step2.set_num_substeps(10);
    
    manager.add_load_step(step1);
    manager.add_load_step(step2);
    
    const LoadStep& retrieved = manager.get_load_step(2);
    EXPECT_EQ(retrieved.id(), 2);
    EXPECT_EQ(retrieved.num_substeps(), 10);
}

TEST(LoadStepManagerTest, GetLoadStepNotFound) {
    LoadStepManager manager;
    
    LoadStep step(1);
    manager.add_load_step(step);
    
    EXPECT_THROW(manager.get_load_step(999), std::out_of_range);
}

TEST(LoadStepManagerTest, Clear) {
    LoadStepManager manager;
    
    LoadStep step(1);
    manager.add_load_step(step);
    
    EXPECT_EQ(manager.num_load_steps(), 1);
    
    manager.clear();
    EXPECT_EQ(manager.num_load_steps(), 0);
}
