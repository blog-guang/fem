#include <gtest/gtest.h>
#include "element/element_base.h"
#include <cmath>

using namespace fem;

// ── Triangle2D ──
TEST(Triangle2DTest, NumNodes) {
    auto e = create_element(ElementType::Triangle2D);
    EXPECT_EQ(e->num_nodes(), 3);
}

TEST(Triangle2DTest, ShapeFunctions_Nodes) {
    auto e = create_element(ElementType::Triangle2D);
    Real N[3];

    // 节点 0: (0,0) → N0=1, N1=0, N2=0
    e->shape_functions({0.0, 0.0, 0.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-15);
    EXPECT_NEAR(N[1], 0.0, 1e-15);
    EXPECT_NEAR(N[2], 0.0, 1e-15);

    // 节点 1: (1,0)
    e->shape_functions({1.0, 0.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-15);
    EXPECT_NEAR(N[1], 1.0, 1e-15);
    EXPECT_NEAR(N[2], 0.0, 1e-15);

    // 节点 2: (0,1)
    e->shape_functions({0.0, 1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-15);
    EXPECT_NEAR(N[1], 0.0, 1e-15);
    EXPECT_NEAR(N[2], 1.0, 1e-15);
}

TEST(Triangle2DTest, ShapeFunctions_PartitionOfUnity) {
    // 任意点形函数之和 = 1
    auto e = create_element(ElementType::Triangle2D);
    Real N[3];

    e->shape_functions({0.2, 0.3, 0.0}, N);
    EXPECT_NEAR(N[0] + N[1] + N[2], 1.0, 1e-15);

    e->shape_functions({0.1, 0.1, 0.0}, N);
    EXPECT_NEAR(N[0] + N[1] + N[2], 1.0, 1e-15);
}

TEST(Triangle2DTest, Gradients_Constant) {
    // 线性三角形梯度为常数
    auto e = create_element(ElementType::Triangle2D);
    Real g1[9], g2[9];

    e->shape_gradients({0.0, 0.0, 0.0}, g1);
    e->shape_gradients({0.3, 0.2, 0.0}, g2);

    for (int i = 0; i < 9; ++i) {
        EXPECT_NEAR(g1[i], g2[i], 1e-15);
    }
}

TEST(Triangle2DTest, QuadPoints) {
    auto e = create_element(ElementType::Triangle2D);
    auto qp = e->quad_points();
    EXPECT_EQ(qp.size, 1);
    EXPECT_NEAR(qp[0].weight, 0.5, 1e-15);  // 参考三角形面积
}

// ── Quad2D ──
TEST(Quad2DTest, NumNodes) {
    auto e = create_element(ElementType::Quad2D);
    EXPECT_EQ(e->num_nodes(), 4);
}

TEST(Quad2DTest, ShapeFunctions_Nodes) {
    auto e = create_element(ElementType::Quad2D);
    Real N[4];

    // 节点 0: (-1,-1)
    e->shape_functions({-1.0, -1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 1.0, 1e-15);
    EXPECT_NEAR(N[1], 0.0, 1e-15);
    EXPECT_NEAR(N[2], 0.0, 1e-15);
    EXPECT_NEAR(N[3], 0.0, 1e-15);

    // 节点 2: (1,1)
    e->shape_functions({1.0, 1.0, 0.0}, N);
    EXPECT_NEAR(N[0], 0.0, 1e-15);
    EXPECT_NEAR(N[1], 0.0, 1e-15);
    EXPECT_NEAR(N[2], 1.0, 1e-15);
    EXPECT_NEAR(N[3], 0.0, 1e-15);
}

TEST(Quad2DTest, ShapeFunctions_PartitionOfUnity) {
    auto e = create_element(ElementType::Quad2D);
    Real N[4];

    e->shape_functions({0.3, -0.5, 0.0}, N);
    EXPECT_NEAR(N[0]+N[1]+N[2]+N[3], 1.0, 1e-15);
}

TEST(Quad2DTest, QuadPoints_WeightSum) {
    // 2x2 Gauss, 权重和 = 4 ([-1,1]x[-1,1] 面积)
    auto e = create_element(ElementType::Quad2D);
    auto qp = e->quad_points();
    EXPECT_EQ(qp.size, 4);

    Real wsum = 0.0;
    for (const auto& q : qp) wsum += q.weight;
    EXPECT_NEAR(wsum, 4.0, 1e-15);
}

// ── Factory ──
TEST(ElementFactoryTest, UnsupportedType) {
    EXPECT_THROW(
        { auto e = create_element(ElementType::Tetrahedron3D); (void)e; },
        std::invalid_argument
    );
}
