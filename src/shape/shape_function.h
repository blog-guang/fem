#pragma once

#include "core/types.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include <memory>
#include <vector>

namespace fem {
namespace shape {

/**
 * ShapeFunction - 形函数抽象基类
 * 
 * 提供统一接口计算：
 * - 形函数值 N(ξ, η, ζ)
 * - 形函数导数 ∂N/∂ξ
 * - 高斯积分点和权重
 * - 雅可比矩阵
 */
class ShapeFunction {
public:
    virtual ~ShapeFunction() = default;
    
    // ═══ 纯虚函数（子类必须实现） ═══
    
    /**
     * 计算形函数值
     * 
     * @param xi 自然坐标 (ξ, η, ζ)
     * @param N  输出：形函数值数组 [N1, N2, ..., Nn]
     */
    virtual void evaluate(const Vec3& xi, Vector& N) const = 0;
    
    /**
     * 计算形函数导数（自然坐标系）
     * 
     * @param xi  自然坐标
     * @param dN  输出：导数矩阵 (nnodes × dim)
     *            dN(i,j) = ∂N_i / ∂ξ_j
     */
    virtual void evaluateDerivatives(const Vec3& xi, DenseMatrix& dN) const = 0;
    
    /**
     * 获取高斯积分点和权重
     * 
     * @param order   积分阶数（1, 2, 3, ...）
     * @param points  输出：积分点坐标
     * @param weights 输出：积分权重
     */
    virtual void getGaussPoints(int order,
                               std::vector<Vec3>& points,
                               std::vector<Real>& weights) const = 0;
    
    /**
     * 获取单元信息
     */
    virtual int dimension() const = 0;
    virtual int numNodes() const = 0;
    virtual ElementType elementType() const = 0;
    
    // ═══ 通用工具函数 ═══
    
    /**
     * 计算雅可比矩阵 J = ∂x/∂ξ
     * 
     * @param xi          自然坐标
     * @param node_coords 节点物理坐标
     * @return 雅可比矩阵 (dim × dim)
     */
    DenseMatrix computeJacobian(const Vec3& xi,
                               const std::vector<Vec3>& node_coords) const;
    
    /**
     * 计算物理坐标系形函数导数
     * 
     * dN/dx = J^{-1} * dN/dξ
     * 
     * @param xi          自然坐标
     * @param node_coords 节点物理坐标
     * @param dN_dx       输出：物理坐标系导数 (nnodes × dim)
     */
    void computePhysicalDerivatives(const Vec3& xi,
                                   const std::vector<Vec3>& node_coords,
                                   DenseMatrix& dN_dx) const;
    
    /**
     * 计算物理坐标（等参插值）
     * 
     * x = Σ N_i * x_i
     * 
     * @param xi          自然坐标
     * @param node_coords 节点物理坐标
     * @return 物理坐标
     */
    Vec3 interpolate(const Vec3& xi,
                    const std::vector<Vec3>& node_coords) const;
};

/**
 * ShapeFunction2D - 2D形函数基类
 */
class ShapeFunction2D : public ShapeFunction {
public:
    int dimension() const override { return 2; }
    
protected:
    // 2D Gauss-Legendre积分点（正方形域）
    void gaussLegendre2D(int order,
                        std::vector<Vec3>& points,
                        std::vector<Real>& weights) const;
    
    // 三角形域高斯积分点
    void gaussTriangle(int order,
                      std::vector<Vec3>& points,
                      std::vector<Real>& weights) const;
};

/**
 * ShapeFunction3D - 3D形函数基类
 */
class ShapeFunction3D : public ShapeFunction {
public:
    int dimension() const override { return 3; }
    
protected:
    // 3D Gauss-Legendre积分点（立方体域）
    void gaussLegendre3D(int order,
                        std::vector<Vec3>& points,
                        std::vector<Real>& weights) const;
    
    // 四面体域高斯积分点
    void gaussTetrahedron(int order,
                         std::vector<Vec3>& points,
                         std::vector<Real>& weights) const;
};

// ═══ 智能指针别名 ═══
using ShapeFunctionPtr = std::shared_ptr<ShapeFunction>;

}  // namespace shape
}  // namespace fem
