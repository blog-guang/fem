/**
 * post_processor_incremental.h - 增量式后处理器
 * 
 * 特点：
 * 1. 增量应力更新（而非重新计算）
 * 2. 状态变量管理（塑性应变、硬化参数等）
 * 3. 真实的本构积分（调用 material->computeStress）
 * 4. 高斯点级别的准确应力
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
#include <unordered_map>

namespace fem {
namespace postprocess {

/**
 * 高斯点状态
 * 存储每个高斯点的应力、应变和状态变量
 */
struct GaussPointState {
    Vector stress;         // 当前应力（Voigt 记号）
    Vector strain;         // 当前总应变
    constitutive::StateVariables state;  // 材料状态变量
    
    GaussPointState(int dimension = 3) {
        int n_comp = (dimension == 3) ? 6 : 3;
        stress.resize(n_comp, 0.0);
        strain.resize(n_comp, 0.0);
    }
};

/**
 * 增量式后处理器
 * 
 * 支持增量分析，维护每个高斯点的完整历史
 */
class IncrementalPostProcessor {
public:
    /**
     * 构造函数
     */
    IncrementalPostProcessor(const Model& model,
                            constitutive::Material* material,
                            int dimension);
    
    /**
     * 初始化（分配所有高斯点的状态）
     */
    void initialize();
    
    /**
     * 增量更新应力和应变
     * 
     * @param displacement_increment 位移增量（当前步的位移）
     * @param data_manager 数据管理器
     * 
     * 功能：
     * 1. 从位移增量计算应变增量
     * 2. 调用 material->computeStress 更新应力
     * 3. 更新状态变量（塑性应变等）
     * 4. 保存到 DataManager
     */
    void update_stress_strain(const Vector& displacement,
                             data::DataManager& data_manager);
    
    /**
     * 重置所有高斯点状态（开始新分析）
     */
    void reset();
    
    /**
     * 提取高斯点应力到 DataManager
     */
    void extract_stress_to_manager(data::DataManager& data_manager,
                                   const std::string& field_name = "stress_gp");
    
    /**
     * 提取高斯点应变到 DataManager
     */
    void extract_strain_to_manager(data::DataManager& data_manager,
                                   const std::string& field_name = "strain_gp");
    
    /**
     * 提取等效塑性应变
     */
    void extract_plastic_strain(data::DataManager& data_manager,
                               const std::string& field_name = "eps_p_eq");
    
    /**
     * 计算 von Mises 应力（从高斯点）
     */
    void compute_von_mises(data::DataManager& data_manager,
                          const std::string& output_field = "von_mises_gp");
    
    /**
     * 外插高斯点数据到节点（加权平均）
     */
    void extrapolate_to_nodes(data::DataManager& data_manager,
                             const std::string& gp_field,
                             const std::string& node_field);
    
    /**
     * 获取高斯点状态（用于调试）
     */
    const GaussPointState& get_gp_state(Index elem_id, Index gp_id) const;
    
private:
    const Model& model_;
    constitutive::Material* material_;
    int dimension_;
    int dofs_per_node_;
    
    // 高斯点状态存储：[elem_id][gp_id] -> GaussPointState
    std::unordered_map<Index, std::vector<GaussPointState>> gp_states_;
    
    // 上一步的位移（用于计算增量）
    Vector displacement_prev_;
    
    /**
     * 计算单元的 B 矩阵（应变-位移关系）
     */
    DenseMatrix compute_B_matrix(const DenseMatrix& dN_dxyz) const;
    
    /**
     * 从单元位移计算单元应变（在高斯点）
     */
    void compute_element_strain_at_gp(Index elem_id,
                                     const Mesh& mesh,
                                     const Vector& displacement,
                                     std::vector<Vector>& strains);
    
    /**
     * 计算 von Mises 应力（从应力向量）
     */
    Real compute_von_mises_stress(const Vector& stress) const;
};

} // namespace postprocess
} // namespace fem
