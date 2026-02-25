#include <gtest/gtest.h>
#include "shape/shape_function_factory.h"
#include <cmath>

using namespace fem;
using namespace fem::shape;

// ═══════════════════════════════════════════════════════════════
// Tri3 测试
// ═══════════════════════════════════════════════════════════════

TEST(Tri3ShapeFunctionTest, NodeValues) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tri3);
    Vector N;
    
    // 节点0: (0, 0)
    shape_func->evaluate(Vec3{0.0, 0.0, 0.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-10);
    EXPECT_NEAR(N[1], 0.0, 1e-10);
    EXPECT_NEAR(N[2], 0.0, 1e-10);
    
    // 节点1: (1, 0)
    shape_func->evaluate(Vec3{1.0, 0.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-10);
    EXPECT_NEAR(N[1], 1.0, 1e-10);
    EXPECT_NEAR(N[2], 0.0, 1e-10);
    
    // 节点2: (0, 1)
    shape_func->evaluate(Vec3{0.0, 1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-10);
    EXPECT_NEAR(N[1], 0.0, 1e-10);
    EXPECT_NEAR(N[2], 1.0, 1e-10);
}

TEST(Tri3ShapeFunctionTest, PartitionOfUnity) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tri3);
    Vector N;
    
    // 中心点
    shape_func->evaluate(Vec3{1.0/3.0, 1.0/3.0, 0.0}, N);
    Real sum = N[0] + N[1] + N[2];
    EXPECT_NEAR(sum, 1.0, 1e-10);
    
    // 任意内部点
    shape_func->evaluate(Vec3{0.2, 0.3, 0.0}, N);
    sum = N[0] + N[1] + N[2];
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

TEST(Tri3ShapeFunctionTest, Derivatives) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tri3);
    DenseMatrix dN;
    
    // 导数是常数（线性单元）
    shape_func->evaluateDerivatives(Vec3{0.5, 0.25, 0.0}, dN);
    
    EXPECT_EQ(dN.rows(), 3);
    EXPECT_EQ(dN.cols(), 2);
    
    // dN1/dξ = -1, dN1/dη = -1
    EXPECT_NEAR(dN(0, 0), -1.0, 1e-10);
    EXPECT_NEAR(dN(0, 1), -1.0, 1e-10);
    
    // dN2/dξ = 1, dN2/dη = 0
    EXPECT_NEAR(dN(1, 0), 1.0, 1e-10);
    EXPECT_NEAR(dN(1, 1), 0.0, 1e-10);
}

TEST(Tri3ShapeFunctionTest, GaussPoints) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tri3);
    std::vector<Vec3> points;
    std::vector<Real> weights;
    
    // 1点积分
    shape_func->getGaussPoints(1, points, weights);
    EXPECT_EQ(points.size(), 1);
    EXPECT_EQ(weights.size(), 1);
    
    // 权重和应为参考三角形面积（0.5）
    Real sum_w = 0.0;
    for (auto w : weights) sum_w += w;
    EXPECT_NEAR(sum_w, 0.5, 1e-10);
}

// ═══════════════════════════════════════════════════════════════
// Quad4 测试
// ═══════════════════════════════════════════════════════════════

TEST(Quad4ShapeFunctionTest, NodeValues) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Quad4);
    Vector N;
    
    // 节点0: (-1, -1)
    shape_func->evaluate(Vec3{-1.0, -1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-10);
    EXPECT_NEAR(N[1], 0.0, 1e-10);
    EXPECT_NEAR(N[2], 0.0, 1e-10);
    EXPECT_NEAR(N[3], 0.0, 1e-10);
    
    // 节点2: (1, 1)
    shape_func->evaluate(Vec3{1.0, 1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-10);
    EXPECT_NEAR(N[1], 0.0, 1e-10);
    EXPECT_NEAR(N[2], 1.0, 1e-10);
    EXPECT_NEAR(N[3], 0.0, 1e-10);
}

TEST(Quad4ShapeFunctionTest, CenterPoint) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Quad4);
    Vector N;
    
    // 中心点 (0, 0)
    shape_func->evaluate(Vec3{0.0, 0.0, 0.0}, N);
    
    // 所有形函数值应相等（对称性）
    EXPECT_NEAR(N[0], 0.25, 1e-10);
    EXPECT_NEAR(N[1], 0.25, 1e-10);
    EXPECT_NEAR(N[2], 0.25, 1e-10);
    EXPECT_NEAR(N[3], 0.25, 1e-10);
}

TEST(Quad4ShapeFunctionTest, GaussPoints) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Quad4);
    std::vector<Vec3> points;
    std::vector<Real> weights;
    
    // 2x2积分
    shape_func->getGaussPoints(2, points, weights);
    EXPECT_EQ(points.size(), 4);
    
    // 权重和应为参考正方形面积（4）
    Real sum_w = 0.0;
    for (auto w : weights) sum_w += w;
    EXPECT_NEAR(sum_w, 4.0, 1e-10);
}

// ═══════════════════════════════════════════════════════════════
// Tet4 测试
// ═══════════════════════════════════════════════════════════════

TEST(Tet4ShapeFunctionTest, NodeValues) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tet4);
    Vector N;
    
    // 节点0: (0, 0, 0)
    shape_func->evaluate(Vec3{0.0, 0.0, 0.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-10);
    EXPECT_NEAR(N[1], 0.0, 1e-10);
    EXPECT_NEAR(N[2], 0.0, 1e-10);
    EXPECT_NEAR(N[3], 0.0, 1e-10);
}

TEST(Tet4ShapeFunctionTest, PartitionOfUnity) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Tet4);
    Vector N;
    
    // 中心点
    shape_func->evaluate(Vec3{0.25, 0.25, 0.25}, N);
    Real sum = N[0] + N[1] + N[2] + N[3];
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// ═══════════════════════════════════════════════════════════════
// Brick8 测试
// ═══════════════════════════════════════════════════════════════

TEST(Brick8ShapeFunctionTest, NodeValues) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Brick8);
    Vector N;
    
    // 节点0: (-1, -1, -1)
    shape_func->evaluate(Vec3{-1.0, -1.0, -1.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-10);
    for (int i = 1; i < 8; ++i) {
        EXPECT_NEAR(N[i], 0.0, 1e-10);
    }
}

TEST(Brick8ShapeFunctionTest, CenterPoint) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Brick8);
    Vector N;
    
    // 中心点 (0, 0, 0)
    shape_func->evaluate(Vec3{0.0, 0.0, 0.0}, N);
    
    // 所有形函数值应相等
    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(N[i], 0.125, 1e-10);
    }
}

// ═══════════════════════════════════════════════════════════════
// 雅可比矩阵测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeFunctionTest, JacobianQuad4) {
    auto shape_func = ShapeFunctionFactory::create(ElementType::Quad4);
    
    // 单位正方形：(0,0), (1,0), (1,1), (0,1)
    std::vector<Vec3> coords = {
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{1.0, 1.0, 0.0},
        Vec3{0.0, 1.0, 0.0}
    };
    
    // 中心点雅可比矩阵
    DenseMatrix J = shape_func->computeJacobian(Vec3{0.0, 0.0, 0.0}, coords);
    
    EXPECT_EQ(J.rows(), 2);
    EXPECT_EQ(J.cols(), 2);
    
    // 单位正方形：J = diag(0.5, 0.5)
    EXPECT_NEAR(J(0, 0), 0.5, 1e-10);
    EXPECT_NEAR(J(1, 1), 0.5, 1e-10);
    EXPECT_NEAR(J(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(J(1, 0), 0.0, 1e-10);
}

// ═══════════════════════════════════════════════════════════════
// 工厂类测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeFunctionFactoryTest, CreateAll) {
    // 测试所有支持的单元类型
    EXPECT_NO_THROW(ShapeFunctionFactory::create(ElementType::Tri3));
    EXPECT_NO_THROW(ShapeFunctionFactory::create(ElementType::Quad4));
    EXPECT_NO_THROW(ShapeFunctionFactory::create(ElementType::Tet4));
    EXPECT_NO_THROW(ShapeFunctionFactory::create(ElementType::Brick8));
}

TEST(ShapeFunctionFactoryTest, ConvenienceAPI) {
    Vector N;
    DenseMatrix dN;
    std::vector<Vec3> points;
    std::vector<Real> weights;
    
    // 便捷接口测试
    ShapeFunctionFactory::evaluate(ElementType::Tri3, Vec3{0.5, 0.25, 0.0}, N);
    EXPECT_EQ(N.size(), 3);
    
    ShapeFunctionFactory::evaluateDerivatives(ElementType::Quad4, Vec3{0.0, 0.0, 0.0}, dN);
    EXPECT_EQ(dN.rows(), 4);
    EXPECT_EQ(dN.cols(), 2);
    
    ShapeFunctionFactory::getGaussPoints(ElementType::Tet4, 1, points, weights);
    EXPECT_EQ(points.size(), 1);
}

// main() provided by gtest_main
