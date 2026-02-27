/**
 * post_processor.h - 后处理计算模块
 * 
 * 功能：
 * - 从位移计算应变
 * - 从应变计算应力
 * - 高斯点数据外插到节点
 * - 高斯点数据平均到单元
 * - 保存到场数据管理器
 */

#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "mesh/model.h"
#include "math/vector.h"
#include "math/dense_matrix.h"
#include "data/data_manager.h"
#include "material/material.h"
#include "shape/shape_function_factory.h"
#include <vector>
#include <memory>

namespace fem {
namespace postprocess {

/**
 * PostProcessor: 后处理计算器
 * 
 * 主要功能：
 * 1. 从位移场计算应变场
 * 2. 从应变场计算应力场
 * 3. 高斯点数据外插到节点
 * 4. 高斯点数据平均到单元中心
 * 5. 保存结果到 DataManager
 * 
 * 使用示例：
 * ```cpp
 * PostProcessor post(model, material, 2);  // 2D
 * 
 * // 计算应变和应力（保存到高斯点）
 * post.compute_strain_stress(displacement, data_manager);
 * 
 * // 外插到节点
 * post.extrapolate_to_nodes(data_manager, "stress_gp", "stress_node");
 * 
 * // 平均到单元
 * post.average_to_elements(data_manager, "stress_gp", "stress_elem");
 * ```
 */
class PostProcessor {
public:
    /**
     * 构造函数
     * @param model 模型引用
     * @param material 材料本构模型
     * @param dimension 维度（2 或 3）
     */
    PostProcessor(const Model& model, 
                  constitutive::Material* material,
                  int dimension);
    
    // ═══════════════════════════════════════════════════════════
    // 主要计算函数
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 计算应变和应力（在高斯点）
     * 
     * @param displacement 位移向量（全局DOF）
     * @param data_manager 数据管理器（输出）
     * @param strain_field_name 应变场名称（默认 "strain_gp"）
     * @param stress_field_name 应力场名称（默认 "stress_gp"）
     * 
     * 输出：
     * - strain_gp: VectorData (6 分量 Voigt 记号)
     * - stress_gp: VectorData (6 分量 Voigt 记号)
     */
    void compute_strain_stress(const Vector& displacement,
                               data::DataManager& data_manager,
                               const std::string& strain_field_name = "strain_gp",
                               const std::string& stress_field_name = "stress_gp");
    
    /**
     * 仅计算应变（在高斯点）
     */
    void compute_strain(const Vector& displacement,
                       data::DataManager& data_manager,
                       const std::string& field_name = "strain_gp");
    
    /**
     * 从应变计算应力（在高斯点）
     */
    void compute_stress_from_strain(data::DataManager& data_manager,
                                    const std::string& strain_field = "strain_gp",
                                    const std::string& stress_field = "stress_gp");
    
    // ═══════════════════════════════════════════════════════════
    // 数据转换函数
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 高斯点数据外插到节点
     * 
     * 使用形函数加权平均：
     * u_node = Σ (N_i * u_gp_i) / Σ N_i
     * 
     * @param data_manager 数据管理器
     * @param gp_field 高斯点场名称
     * @param node_field 节点场名称（输出）
     */
    void extrapolate_to_nodes(data::DataManager& data_manager,
                             const std::string& gp_field,
                             const std::string& node_field);
    
    /**
     * 高斯点数据平均到单元中心
     * 
     * u_elem = (1/n_gp) * Σ u_gp_i
     * 
     * @param data_manager 数据管理器
     * @param gp_field 高斯点场名称
     * @param elem_field 单元场名称（输出）
     */
    void average_to_elements(data::DataManager& data_manager,
                            const std::string& gp_field,
                            const std::string& elem_field);
    
    // ═══════════════════════════════════════════════════════════
    // 应力分量提取（用于可视化）
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 提取应力分量到标量场
     * 
     * @param data_manager 数据管理器
     * @param stress_field 应力场名称（VectorData，6 分量）
     * @param component 分量索引（0=σ_xx, 1=σ_yy, 2=σ_zz, 3=τ_xy, 4=τ_yz, 5=τ_xz）
     * @param output_field 输出标量场名称
     */
    void extract_stress_component(data::DataManager& data_manager,
                                  const std::string& stress_field,
                                  int component,
                                  const std::string& output_field);
    
    /**
     * 计算 von Mises 等效应力
     * 
     * σ_vm = sqrt(1.5 * s:s)
     * s = σ - (1/3)*tr(σ)*I  (偏应力)
     * 
     * @param data_manager 数据管理器
     * @param stress_field 应力场名称
     * @param output_field 输出标量场名称（默认 "von_mises"）
     */
    void compute_von_mises(data::DataManager& data_manager,
                          const std::string& stress_field,
                          const std::string& output_field = "von_mises");
    
    // ═══════════════════════════════════════════════════════════
    // 工具函数
    // ═══════════════════════════════════════════════════════════
    
    /**
     * 获取每单元高斯点数量
     */
    Index get_num_gauss_points(ElementType elem_type) const;
    
    /**
     * 获取应变/应力分量数量（2D=3, 3D=6）
     */
    Index num_strain_components() const { return dimension_ == 2 ? 3 : 6; }

private:
    const Model& model_;                  ///< 模型引用
    constitutive::Material* material_;    ///< 材料本构
    int dimension_;                       ///< 维度（2 或 3）
    Index dofs_per_node_;                 ///< 每节点自由度数
    
    /**
     * 计算单元的应变（在高斯点）
     */
    void compute_element_strain(Index elem_id,
                               const Mesh& mesh,
                               const Vector& displacement,
                               std::vector<Vector>& strains);
    
    /**
     * 从应变计算应力（单个高斯点）
     */
    void compute_stress_at_gp(const Vector& strain,
                             Vector& stress,
                             constitutive::StateVariables& state);
    
    /**
     * 构造 B 矩阵（2D）
     */
    DenseMatrix build_B_matrix_2D(const DenseMatrix& dN_dxyz) const;
    
    /**
     * 构造 B 矩阵（3D）
     */
    DenseMatrix build_B_matrix_3D(const DenseMatrix& dN_dxyz) const;
    
    /**
     * 计算单元位移向量
     */
    Vector get_element_displacement(Index elem_id,
                                   const Mesh& mesh,
                                   const Vector& displacement) const;
};

} // namespace postprocess
} // namespace fem
