#pragma once

#include "core/types.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "material/state_variables.h"
#include <string>
#include <map>
#include <memory>

namespace fem {
namespace constitutive {

/**
 * Material: 材料本构模型基类
 * 
 * 定义材料本构模型的统一接口，子类实现具体材料模型：
 * - 各向同性/各向异性弹性
 * - 塑性（J2, Drucker-Prager等）
 * - 粘塑性
 * - 超弹性
 * - 损伤材料
 * 
 * 核心职责：
 * 1. 应力更新（给定应变增量，计算应力）
 * 2. 切线刚度矩阵计算（D_mat: dσ/dε）
 * 3. 几何刚度矩阵计算（K_geo，用于非线性分析）
 * 4. 状态变量管理（塑性应变、损伤等）
 */
class Material {
public:
    // ═══ 构造与析构 ═══
    virtual ~Material() = default;
    
    // ═══ 核心接口 (子类必须实现) ═══
    
    /**
     * 应力更新：给定应变增量，计算新的应力
     * 
     * @param strain_inc  应变增量 (Voigt记号: [ε11, ε22, ε33, γ12, γ23, γ13])
     * @param stress      输出: 更新后的应力 (Voigt记号)
     * @param state       输入/输出: 材料状态变量（会被修改）
     * 
     * 备注：
     * - 对于弹性材料，state不变
     * - 对于塑性材料，需更新塑性应变、硬化变量等
     * - 应变使用小应变假设（几何线性）或有限应变（几何非线性）
     */
    virtual void computeStress(
        const Vector& strain_inc, 
        Vector& stress, 
        StateVariables& state
    ) = 0;
    
    /**
     * 计算切线刚度矩阵（材料刚度）
     * 
     * @param strain  当前总应变
     * @param D_mat   输出: 切线刚度矩阵 (size: n×n, Voigt记号)
     * @param state   当前状态变量（只读）
     * 
     * 备注：
     * - 弹性材料：D_mat = 弹性刚度矩阵（常数）
     * - 塑性材料：D_mat = 弹塑性切线刚度（与应力状态有关）
     * - 几何非线性：可能需要当前应力信息
     */
    virtual void computeTangent(
        const Vector& strain,
        DenseMatrix& D_mat,
        const StateVariables& state
    ) = 0;
    
    /**
     * 计算几何刚度矩阵（用于几何非线性分析）
     * 
     * @param stress  当前应力
     * @param K_geo   输出: 几何刚度矩阵
     * 
     * 备注：
     * - 线性分析中通常为零
     * - 大变形/屈曲分析中必须考虑
     * - 具体形式与单元类型和应力状态相关
     */
    virtual void computeGeometricStiffness(
        const Vector& stress,
        DenseMatrix& K_geo
    );
    
    /**
     * 计算应变能密度
     * 
     * @param strain  当前应变
     * @param state   当前状态变量
     * @return 应变能密度 Ψ(ε)
     */
    virtual Real strainEnergy(
        const Vector& strain,
        const StateVariables& state
    ) const = 0;
    
    // ═══ 状态管理 ═══
    
    /**
     * 创建初始状态变量（工厂方法）
     * 
     * @return 初始化的状态变量对象
     */
    virtual StateVariables createState() const = 0;
    
    /**
     * 初始化状态变量（为已有对象设置初值）
     * 
     * @param state  要初始化的状态变量
     */
    virtual void initializeState(StateVariables& state) const;
    
    // ═══ 参数管理 ═══
    
    /**
     * 设置材料参数
     * 
     * @param name  参数名 (如 "E", "nu", "sigma_y")
     * @param value 参数值
     * @throws std::invalid_argument 如果参数名不存在或值无效
     */
    virtual void setParameter(const std::string& name, Real value);
    
    /**
     * 获取材料参数
     * 
     * @param name  参数名
     * @return 参数值
     * @throws std::invalid_argument 如果参数名不存在
     */
    virtual Real getParameter(const std::string& name) const;
    
    /**
     * 获取所有参数名列表
     */
    virtual std::vector<std::string> getParameterNames() const;
    
    // ═══ 工具接口 ═══
    
    /**
     * 获取材料类型名称（用于调试/日志）
     */
    virtual std::string typeName() const = 0;
    
    /**
     * 获取应变/应力分量数（Voigt记号）
     * - 2D平面应力/应变: 3 (ε11, ε22, γ12)
     * - 2D平面应变/轴对称: 4 (ε11, ε22, ε33, γ12)
     * - 3D: 6 (ε11, ε22, ε33, γ12, γ23, γ13)
     */
    std::size_t strainSize() const { return strain_size_; }
    
    /**
     * 打印材料信息
     */
    virtual void print() const;
    
    /**
     * 验证参数合法性（子类可重写）
     * 
     * @throws std::invalid_argument 如果参数非法
     */
    virtual void validateParameters() const;
    
protected:
    // ═══ 子类访问成员 ═══
    std::size_t strain_size_{6};  // 默认3D (6分量)
    std::map<std::string, Real> parameters_;  // 材料参数存储
    
    // 构造函数（仅供子类）
    explicit Material(std::size_t strain_size = 6) 
        : strain_size_(strain_size) {}
};

// ═══ 智能指针别名 ═══
using MaterialPtr = std::shared_ptr<Material>;

}  // namespace constitutive
}  // namespace fem
