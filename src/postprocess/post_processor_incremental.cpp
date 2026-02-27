/**
 * post_processor_incremental.cpp - 增量式后处理器实现
 */

#include "postprocess/post_processor_incremental.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace postprocess {

IncrementalPostProcessor::IncrementalPostProcessor(const Model& model,
                                                 constitutive::Material* material,
                                                 int dimension)
    : model_(model), material_(material), dimension_(dimension) {
    
    dofs_per_node_ = dimension;
}

void IncrementalPostProcessor::initialize() {
    gp_states_.clear();
    
    // 为每个网格的每个单元的每个高斯点分配状态
    for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            const Element& elem = mesh.element(elem_id);
            
            // 创建形函数对象以获取高斯点数量
            auto shape = shape::ShapeFunctionFactory::create(elem.type());
            
            std::vector<Vec3> gp_coords;
            std::vector<Real> gp_weights;
            shape->getGaussPoints(2, gp_coords, gp_weights);
            
            std::vector<GaussPointState> gp_list;
            for (size_t gp = 0; gp < gp_coords.size(); ++gp) {
                GaussPointState gp_state(dimension_);
                
                // 初始化状态变量
                gp_state.state = material_->createState();
                
                gp_list.push_back(gp_state);
            }
            
            gp_states_[elem_id] = gp_list;
        }
    }
    
    FEM_INFO("IncrementalPostProcessor initialized: " + 
            std::to_string(gp_states_.size()) + " elements");
}

void IncrementalPostProcessor::reset() {
    for (auto& [elem_id, gp_list] : gp_states_) {
        for (auto& gp_state : gp_list) {
            // 重置应力、应变
            for (size_t i = 0; i < gp_state.stress.size(); ++i) {
                gp_state.stress[i] = 0.0;
                gp_state.strain[i] = 0.0;
            }
            
            // 重置状态变量
            gp_state.state = material_->createState();
        }
    }
    
    displacement_prev_.resize(0);
}

void IncrementalPostProcessor::update_stress_strain(const Vector& displacement,
                                                   data::DataManager& data_manager) {
    // 如果是第一步，初始化
    if (displacement_prev_.size() == 0) {
        displacement_prev_.resize(displacement.size(), 0.0);
    }
    
    // 计算位移增量
    Vector displacement_increment(displacement.size());
    for (size_t i = 0; i < displacement.size(); ++i) {
        displacement_increment[i] = displacement[i] - displacement_prev_[i];
    }
    
    // 遍历所有网格
    for (std::size_t mesh_id = 0; mesh_id < model_.num_meshes(); ++mesh_id) {
        const Mesh& mesh = model_.mesh(mesh_id);
        
        // 遍历所有单元
        for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            const Element& elem = mesh.element(elem_id);
            
            // 获取单元节点坐标
            std::vector<Vec3> coords;
            for (Index node_id : elem.nodes()) {
                coords.push_back(mesh.node(node_id).coords());
            }
            
            // 创建形函数对象
            auto shape = shape::ShapeFunctionFactory::create(elem.type());
            
            // 获取高斯积分点
            std::vector<Vec3> gp_coords;
            std::vector<Real> gp_weights;
            shape->getGaussPoints(2, gp_coords, gp_weights);  // 2阶积分
            
            // 提取单元位移增量
            Index n_nodes = elem.nodes().size();
            Vector u_elem_inc(n_nodes * dofs_per_node_);
            for (Index i = 0; i < n_nodes; ++i) {
                Index node_id = elem.nodes()[i];
                for (Index d = 0; d < dofs_per_node_; ++d) {
                    u_elem_inc[i * dofs_per_node_ + d] = displacement_increment[node_id * dofs_per_node_ + d];
                }
            }
            
            // 在每个高斯点计算应变增量和更新应力
            for (size_t gp_idx = 0; gp_idx < gp_coords.size(); ++gp_idx) {
                const Vec3& xi = gp_coords[gp_idx];
                
                // 计算形函数导数（物理坐标系）
                DenseMatrix dN_dxyz;
                shape->computePhysicalDerivatives(xi, coords, dN_dxyz);
                
                // 构造 B 矩阵
                DenseMatrix B = compute_B_matrix(dN_dxyz);
                
                // 应变增量 = B * u_increment
                Vector strain_increment(B.rows(), 0.0);
                for (std::size_t i = 0; i < B.rows(); ++i) {
                    for (std::size_t j = 0; j < B.cols(); ++j) {
                        strain_increment[i] += B(i, j) * u_elem_inc[j];
                    }
                }
                
                // 获取高斯点状态
                GaussPointState& gp_state = gp_states_[elem_id][gp_idx];
                
                // 更新总应变
                for (size_t i = 0; i < strain_increment.size(); ++i) {
                    gp_state.strain[i] += strain_increment[i];
                }
                
                // 调用材料本构更新应力
                material_->computeStress(strain_increment, gp_state.stress, gp_state.state);
            }
        }
    }
    
    // 保存当前位移
    displacement_prev_ = displacement;
}

DenseMatrix IncrementalPostProcessor::compute_B_matrix(const DenseMatrix& dN_dxyz) const {
    Index n_nodes = dN_dxyz.cols();
    
    if (dimension_ == 2) {
        // 2D: B 矩阵 (3 × 2*n_nodes)
        DenseMatrix B(3, 2 * n_nodes, 0.0);
        
        for (Index i = 0; i < n_nodes; ++i) {
            Real dN_dx = dN_dxyz(0, i);
            Real dN_dy = dN_dxyz(1, i);
            
            // ε_xx = du/dx
            B(0, 2 * i) = dN_dx;
            
            // ε_yy = dv/dy
            B(1, 2 * i + 1) = dN_dy;
            
            // γ_xy = du/dy + dv/dx
            B(2, 2 * i) = dN_dy;
            B(2, 2 * i + 1) = dN_dx;
        }
        
        return B;
        
    } else {
        // 3D: B 矩阵 (6 × 3*n_nodes)
        DenseMatrix B(6, 3 * n_nodes, 0.0);
        
        for (Index i = 0; i < n_nodes; ++i) {
            Real dN_dx = dN_dxyz(0, i);
            Real dN_dy = dN_dxyz(1, i);
            Real dN_dz = dN_dxyz(2, i);
            
            // ε_xx = du/dx
            B(0, 3 * i) = dN_dx;
            
            // ε_yy = dv/dy
            B(1, 3 * i + 1) = dN_dy;
            
            // ε_zz = dw/dz
            B(2, 3 * i + 2) = dN_dz;
            
            // γ_xy = du/dy + dv/dx
            B(3, 3 * i) = dN_dy;
            B(3, 3 * i + 1) = dN_dx;
            
            // γ_yz = dv/dz + dw/dy
            B(4, 3 * i + 1) = dN_dz;
            B(4, 3 * i + 2) = dN_dy;
            
            // γ_xz = du/dz + dw/dx
            B(5, 3 * i) = dN_dz;
            B(5, 3 * i + 2) = dN_dx;
        }
        
        return B;
    }
}

void IncrementalPostProcessor::extract_stress_to_manager(data::DataManager& data_manager,
                                                        const std::string& field_name) {
    // 统计高斯点总数
    Index total_gp = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        total_gp += gp_list.size();
    }
    
    int n_comp = (dimension_ == 3) ? 6 : 3;
    
    // 创建场数据
    auto* stress_field = data_manager.add_field<data::VectorData>(
        field_name, data::DataLocation::GaussPoint, total_gp,
        Vector(n_comp, 0.0));
    
    // 填充数据
    Index gp_global_id = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        for (const auto& gp_state : gp_list) {
            stress_field->set(gp_global_id, gp_state.stress);
            gp_global_id++;
        }
    }
}

void IncrementalPostProcessor::extract_strain_to_manager(data::DataManager& data_manager,
                                                        const std::string& field_name) {
    Index total_gp = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        total_gp += gp_list.size();
    }
    
    int n_comp = (dimension_ == 3) ? 6 : 3;
    
    auto* strain_field = data_manager.add_field<data::VectorData>(
        field_name, data::DataLocation::GaussPoint, total_gp,
        Vector(n_comp, 0.0));
    
    Index gp_global_id = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        for (const auto& gp_state : gp_list) {
            strain_field->set(gp_global_id, gp_state.strain);
            gp_global_id++;
        }
    }
}

void IncrementalPostProcessor::extract_plastic_strain(data::DataManager& data_manager,
                                                     const std::string& field_name) {
    Index total_gp = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        total_gp += gp_list.size();
    }
    
    auto* eps_p_field = data_manager.add_field<data::RealData>(
        field_name, data::DataLocation::GaussPoint, total_gp, 0.0);
    
    Index gp_global_id = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        for (const auto& gp_state : gp_list) {
            eps_p_field->set(gp_global_id, gp_state.state.equiv_plastic_strain);
            gp_global_id++;
        }
    }
}

Real IncrementalPostProcessor::compute_von_mises_stress(const Vector& stress) const {
    if (dimension_ == 3) {
        // 3D: σ_vm = √(0.5*[(σ_xx-σ_yy)² + (σ_yy-σ_zz)² + (σ_zz-σ_xx)²] + 3*(τ_xy² + τ_yz² + τ_xz²))
        Real s_xx = stress[0];
        Real s_yy = stress[1];
        Real s_zz = stress[2];
        Real s_xy = stress[3];
        Real s_yz = stress[4];
        Real s_xz = stress[5];
        
        Real diff1 = s_xx - s_yy;
        Real diff2 = s_yy - s_zz;
        Real diff3 = s_zz - s_xx;
        
        Real vm_sq = 0.5 * (diff1*diff1 + diff2*diff2 + diff3*diff3) +
                    3.0 * (s_xy*s_xy + s_yz*s_yz + s_xz*s_xz);
        
        return std::sqrt(vm_sq);
        
    } else {
        // 2D: σ_vm = √(σ_xx² + σ_yy² - σ_xx*σ_yy + 3*τ_xy²)
        Real s_xx = stress[0];
        Real s_yy = stress[1];
        Real s_xy = stress[2];
        
        Real vm_sq = s_xx*s_xx + s_yy*s_yy - s_xx*s_yy + 3.0*s_xy*s_xy;
        
        return std::sqrt(vm_sq);
    }
}

void IncrementalPostProcessor::compute_von_mises(data::DataManager& data_manager,
                                                const std::string& output_field) {
    Index total_gp = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        total_gp += gp_list.size();
    }
    
    auto* vm_field = data_manager.add_field<data::RealData>(
        output_field, data::DataLocation::GaussPoint, total_gp, 0.0);
    
    Index gp_global_id = 0;
    for (const auto& [elem_id, gp_list] : gp_states_) {
        for (const auto& gp_state : gp_list) {
            Real vm = compute_von_mises_stress(gp_state.stress);
            vm_field->set(gp_global_id, vm);
            gp_global_id++;
        }
    }
}

const GaussPointState& IncrementalPostProcessor::get_gp_state(Index elem_id, Index gp_id) const {
    return gp_states_.at(elem_id).at(gp_id);
}

} // namespace postprocess
} // namespace fem
