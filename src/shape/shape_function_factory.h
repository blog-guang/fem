#pragma once

#include "shape_function.h"
#include "shape_tri3.h"
#include "shape_quad4.h"
#include "shape_tet4.h"
#include "shape_brick8.h"
#include <memory>

namespace fem {
namespace shape {

/**
 * ShapeFunctionFactory - 形函数工厂类
 * 
 * 根据单元类型创建对应的形函数对象
 */
class ShapeFunctionFactory {
public:
    /**
     * 创建形函数对象
     * 
     * @param type 单元类型
     * @return 形函数智能指针
     */
    static ShapeFunctionPtr create(ElementType type);
    
    /**
     * 便捷接口：计算形函数值
     * 
     * @param type 单元类型
     * @param xi   自然坐标
     * @param N    输出：形函数值
     */
    static void evaluate(ElementType type, const Vec3& xi, Vector& N);
    
    /**
     * 便捷接口：计算形函数导数
     * 
     * @param type 单元类型
     * @param xi   自然坐标
     * @param dN   输出：形函数导数
     */
    static void evaluateDerivatives(ElementType type, const Vec3& xi, DenseMatrix& dN);
    
    /**
     * 便捷接口：获取高斯积分点
     * 
     * @param type    单元类型
     * @param order   积分阶数
     * @param points  输出：积分点
     * @param weights 输出：积分权重
     */
    static void getGaussPoints(ElementType type, int order,
                              std::vector<Vec3>& points,
                              std::vector<Real>& weights);
};

}  // namespace shape
}  // namespace fem
