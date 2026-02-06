#include <gtest/gtest.h>
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "math/sparse_matrix.h"

using namespace fem;

// ═══ Vector Tests ═══
TEST(VectorTest, Construction) {
    Vector v1;
    EXPECT_EQ(v1.size(), 0u);
    
    Vector v2(5, 2.0);
    EXPECT_EQ(v2.size(), 5u);
    EXPECT_DOUBLE_EQ(v2[0], 2.0);
    EXPECT_DOUBLE_EQ(v2[4], 2.0);
}

TEST(VectorTest, Arithmetic) {
    Vector v1(3, 1.0);
    Vector v2(3, 2.0);
    
    Vector v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3[0], 3.0);
    
    Vector v4 = v2 - v1;
    EXPECT_DOUBLE_EQ(v4[0], 1.0);
    
    Vector v5 = v1 * 3.0;
    EXPECT_DOUBLE_EQ(v5[0], 3.0);
    
    Vector v6 = v2 / 2.0;
    EXPECT_DOUBLE_EQ(v6[0], 1.0);
}

TEST(VectorTest, DotProduct) {
    Vector v1(3);
    v1[0] = 1.0; v1[1] = 2.0; v1[2] = 3.0;
    
    Vector v2(3);
    v2[0] = 4.0; v2[1] = 5.0; v2[2] = 6.0;
    
    Real dot = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dot, 1*4 + 2*5 + 3*6);  // 32
}

TEST(VectorTest, Norms) {
    Vector v(3);
    v[0] = 3.0; v[1] = 4.0; v[2] = 0.0;
    
    EXPECT_DOUBLE_EQ(v.norm(), 5.0);         // sqrt(9+16)
    EXPECT_DOUBLE_EQ(v.norm1(), 7.0);        // 3+4
    EXPECT_DOUBLE_EQ(v.norm_inf(), 4.0);     // max(3,4,0)
    EXPECT_DOUBLE_EQ(v.norm_squared(), 25.0);
}

TEST(VectorTest, Normalize) {
    Vector v(2);
    v[0] = 3.0; v[1] = 4.0;
    
    v.normalize();
    EXPECT_NEAR(v.norm(), 1.0, 1e-12);
    EXPECT_NEAR(v[0], 0.6, 1e-12);  // 3/5
    EXPECT_NEAR(v[1], 0.8, 1e-12);  // 4/5
}

// ═══ DenseMatrix Tests ═══
TEST(DenseMatrixTest, Construction) {
    DenseMatrix A;
    EXPECT_EQ(A.rows(), 0u);
    EXPECT_EQ(A.cols(), 0u);
    
    DenseMatrix B(3, 4, 2.0);
    EXPECT_EQ(B.rows(), 3u);
    EXPECT_EQ(B.cols(), 4u);
    EXPECT_DOUBLE_EQ(B(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(B(2, 3), 2.0);
}

TEST(DenseMatrixTest, Identity) {
    DenseMatrix I(3, 3);
    I.identity();
    
    EXPECT_DOUBLE_EQ(I(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(I(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(I(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(I(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(I(1, 2), 0.0);
}

TEST(DenseMatrixTest, Arithmetic) {
    DenseMatrix A(2, 2);
    A(0, 0) = 1; A(0, 1) = 2;
    A(1, 0) = 3; A(1, 1) = 4;
    
    DenseMatrix B(2, 2);
    B(0, 0) = 5; B(0, 1) = 6;
    B(1, 0) = 7; B(1, 1) = 8;
    
    DenseMatrix C = A + B;
    EXPECT_DOUBLE_EQ(C(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(C(1, 1), 12.0);
    
    DenseMatrix D = B - A;
    EXPECT_DOUBLE_EQ(D(0, 0), 4.0);
    
    DenseMatrix E = A * 2.0;
    EXPECT_DOUBLE_EQ(E(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(E(1, 1), 8.0);
}

TEST(DenseMatrixTest, MatrixMultiplication) {
    DenseMatrix A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    
    DenseMatrix B(3, 2);
    B(0, 0) = 7; B(0, 1) = 8;
    B(1, 0) = 9; B(1, 1) = 10;
    B(2, 0) = 11; B(2, 1) = 12;
    
    DenseMatrix C = A * B;  // 2x2
    EXPECT_EQ(C.rows(), 2u);
    EXPECT_EQ(C.cols(), 2u);
    
    // C[0,0] = 1*7 + 2*9 + 3*11 = 58
    EXPECT_DOUBLE_EQ(C(0, 0), 58.0);
    // C[0,1] = 1*8 + 2*10 + 3*12 = 64
    EXPECT_DOUBLE_EQ(C(0, 1), 64.0);
    // C[1,0] = 4*7 + 5*9 + 6*11 = 139
    EXPECT_DOUBLE_EQ(C(1, 0), 139.0);
}

TEST(DenseMatrixTest, MatVec) {
    DenseMatrix A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    
    Vector x(3);
    x[0] = 1; x[1] = 2; x[2] = 3;
    
    Vector y = A * x;
    EXPECT_EQ(y.size(), 2u);
    EXPECT_DOUBLE_EQ(y[0], 1*1 + 2*2 + 3*3);  // 14
    EXPECT_DOUBLE_EQ(y[1], 4*1 + 5*2 + 6*3);  // 32
}

TEST(DenseMatrixTest, Transpose) {
    DenseMatrix A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    
    DenseMatrix B = A.transpose();
    EXPECT_EQ(B.rows(), 3u);
    EXPECT_EQ(B.cols(), 2u);
    EXPECT_DOUBLE_EQ(B(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(B(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(B(2, 1), 6.0);
}

TEST(DenseMatrixTest, Symmetric) {
    DenseMatrix A(3, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 2; A(1, 1) = 4; A(1, 2) = 5;
    A(2, 0) = 3; A(2, 1) = 5; A(2, 2) = 6;
    
    EXPECT_TRUE(A.is_symmetric());
    
    A(0, 1) = 2.1;
    EXPECT_FALSE(A.is_symmetric());
}

// ═══ SparseMatrix Tests ═══
TEST(SparseMatrixTest, COO_Construction) {
    SparseMatrixCOO coo(3, 3);
    coo.add(0, 0, 1.0);
    coo.add(0, 1, 2.0);
    coo.add(1, 1, 3.0);
    coo.add(2, 2, 4.0);
    
    EXPECT_EQ(coo.rows(), 3u);
    EXPECT_EQ(coo.cols(), 3u);
    EXPECT_EQ(coo.nnz(), 4u);
}

TEST(SparseMatrixTest, COO_to_CSR) {
    SparseMatrixCOO coo(3, 3);
    coo.add(0, 0, 1.0);
    coo.add(0, 2, 2.0);
    coo.add(1, 1, 3.0);
    coo.add(2, 0, 4.0);
    coo.add(2, 2, 5.0);
    
    SparseMatrixCSR csr = coo_to_csr(coo);
    EXPECT_EQ(csr.rows(), 3u);
    EXPECT_EQ(csr.cols(), 3u);
    EXPECT_EQ(csr.nnz(), 5u);
    
    // 检查 row_ptr
    EXPECT_EQ(csr.row_ptr()[0], 0u);
    EXPECT_EQ(csr.row_ptr()[1], 2u);  // 第0行有2个非零元
    EXPECT_EQ(csr.row_ptr()[2], 3u);  // 第1行有1个非零元
    EXPECT_EQ(csr.row_ptr()[3], 5u);  // 第2行有2个非零元
}

TEST(SparseMatrixTest, CSR_MatVec) {
    // 创建简单的 2x2 单位矩阵
    SparseMatrixCOO coo(2, 2);
    coo.add(0, 0, 2.0);
    coo.add(0, 1, 1.0);
    coo.add(1, 0, 1.0);
    coo.add(1, 1, 3.0);
    
    SparseMatrixCSR csr = coo_to_csr(coo);
    
    Vector x(2);
    x[0] = 1.0;
    x[1] = 2.0;
    
    Vector y = csr * x;
    EXPECT_EQ(y.size(), 2u);
    EXPECT_DOUBLE_EQ(y[0], 2*1 + 1*2);  // 4
    EXPECT_DOUBLE_EQ(y[1], 1*1 + 3*2);  // 7
}

TEST(SparseMatrixTest, COO_Duplicate_Merge) {
    // 测试重复元素合并
    SparseMatrixCOO coo(2, 2);
    coo.add(0, 0, 1.0);
    coo.add(0, 0, 2.0);  // 重复
    coo.add(1, 1, 3.0);
    
    SparseMatrixCSR csr = coo_to_csr(coo);
    EXPECT_EQ(csr.nnz(), 2u);  // 合并后只有2个非零元
    
    // 第一个元素应该是 1+2=3
    Vector x(2, 1.0);
    Vector y = csr * x;
    EXPECT_DOUBLE_EQ(y[0], 3.0);  // (0,0) = 1+2 = 3
    EXPECT_DOUBLE_EQ(y[1], 3.0);  // (1,1) = 3
}

TEST(SparseMatrixTest, COO_to_CSC) {
    SparseMatrixCOO coo(3, 3);
    coo.add(0, 0, 1.0);
    coo.add(1, 0, 2.0);
    coo.add(2, 1, 3.0);
    coo.add(0, 2, 4.0);
    
    SparseMatrixCSC csc = coo_to_csc(coo);
    EXPECT_EQ(csc.rows(), 3u);
    EXPECT_EQ(csc.cols(), 3u);
    EXPECT_EQ(csc.nnz(), 4u);
    
    // 检查 col_ptr
    EXPECT_EQ(csc.col_ptr()[0], 0u);
    EXPECT_EQ(csc.col_ptr()[1], 2u);  // 第0列有2个非零元
    EXPECT_EQ(csc.col_ptr()[2], 3u);  // 第1列有1个非零元
    EXPECT_EQ(csc.col_ptr()[3], 4u);  // 第2列有1个非零元
}

TEST(SparseMatrixTest, EmptyMatrix) {
    SparseMatrixCOO coo(5, 5);
    // 空矩阵
    
    SparseMatrixCSR csr = coo_to_csr(coo);
    EXPECT_EQ(csr.rows(), 5u);
    EXPECT_EQ(csr.cols(), 5u);
    EXPECT_EQ(csr.nnz(), 0u);
    
    Vector x(5, 1.0);
    Vector y = csr * x;
    EXPECT_EQ(y.size(), 5u);
    for (std::size_t i = 0; i < 5; ++i) {
        EXPECT_DOUBLE_EQ(y[i], 0.0);
    }
}
