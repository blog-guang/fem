#include "postprocess/post_processor.h"
#include "physics/physics_base.h"
#include "core/logger.h"
#include <cmath>

namespace fem {
namespace postprocess {

using namespace physics;

PostProcessor::PostProcessor(const Model& model,
                            constitutive::Material* material,
                            int dimension)
    : model_(model), material_(material), dimension_(dimension) {
    
    dofs_per_node_ = dimension_;
    
    if (!material_) {
        FEM_ERROR("PostProcessor: material is null");
    }
    if (dimension_ != 2 && dimension_ != 3) {
        FEM_ERROR("PostProcessor: dimension must be 2 or 3");
    }
}

// ═══════════════════════════════════════════════════════════
// 主要计算函数
// ═══════════════════════════════════════════════════════════

void PostProcessor::compute_strain_stress(const Vector& displacement,
                                         data::DataManager& data_manager,
                                         const std::string& strain_field_name,
                                         const std::string& stress_field_name) {
    // 计算应变
    compute_strain(displacement, data_manager, strain_field_name);
    
    // 从应变计算应力
    compute_stress_from_strain(data_manager, strain_field_name, stress_field_name);
}

void PostProcessor::compute_strain(const Vector& displacement,
                                  data::DataManager& data_manager,
                                  const std::string& field_name) {
    // 获取第一个网格（假设只有一个网格）
    if (model_.num_meshes() == 0) {
        FEM_WARN("PostProcessor: no mesh in model");
        return;
    }
    
    const Mesh& mesh = model_.mesh(0);
    Index n_elems = mesh.num_elements();
    
    // 计算总高斯点数
    Index total_gp = 0;
    for (Index elem_id = 0; elem_id < n_elems; ++elem_id) {
        ElementType elem_type = mesh.element(elem_id).type();
        Index n_gp = get_num_gauss_points(elem_type);
        total_gp += n_gp;
    }
    
    // 创建应变场（VectorData，每个高斯点存储 3 或 6 个分量）
    Index num_components = num_strain_components();
    
    auto* strain_field = data_manager.add_field<data::VectorData>(
        field_name, data::DataLocation::GaussPoint, total_gp,
        Vector(num_components, 0.0));  // 默认值：零向量
    
    // 遍历所有单元
    Index gp_offset = 0;
    for (Index elem_id = 0; elem_id < n_elems; ++elem_id) {
        // 计算单元的应变（在所有高斯点）
        std::vector<Vector> strains;
        compute_element_strain(elem_id, mesh, displacement, strains);
        
        // 保存到场数据
        for (const auto& strain : strains) {
            strain_field->set(gp_offset++, strain);
        }
    }
    
    FEM_INFO("Computed strain field at " + std::to_string(total_gp) + " Gauss points");
}

void PostProcessor::compute_stress_from_strain(data::DataManager& data_manager,
                                              const std::string& strain_field,
                                              const std::string& stress_field) {
    // 获取应变场
    auto* strain_data = data_manager.get_field<data::VectorData>(strain_field);
    if (!strain_data) {
        FEM_ERROR("PostProcessor: strain field '" + strain_field + "' not found");
        return;
    }
    
    Index total_gp = strain_data->size();
    
    // 创建应力场
    auto* stress_data = data_manager.add_field<data::VectorData>(
        stress_field, data::DataLocation::GaussPoint, total_gp);
    
    // 创建状态变量
    constitutive::StateVariables state = material_->createState();
    
    // 遍历所有高斯点
    for (Index gp = 0; gp < total_gp; ++gp) {
        const Vector& strain = strain_data->get(gp);
        
        Vector stress(num_strain_components());
        compute_stress_at_gp(strain, stress, state);
        
        stress_data->set(gp, stress);
    }
    
    FEM_INFO("Computed stress field at " + std::to_string(total_gp) + " Gauss points");
}

// ═══════════════════════════════════════════════════════════
// 数据转换函数
// ═══════════════════════════════════════════════════════════

void PostProcessor::extrapolate_to_nodes(data::DataManager& data_manager,
                                        const std::string& gp_field,
                                        const std::string& node_field) {
    // 获取高斯点场
    auto* gp_data = data_manager.get_field<data::VectorData>(gp_field);
    if (!gp_data) {
        FEM_ERROR("PostProcessor: Gauss point field '" + gp_field + "' not found");
        return;
    }
    
    const Mesh& mesh = model_.mesh(0);
    Index n_nodes = mesh.num_nodes();
    Index n_elems = mesh.num_elements();
    
    // 创建节点场
    auto* node_data = data_manager.add_field<data::VectorData>(
        node_field, data::DataLocation::Node, n_nodes);
    
    // 累加器（用于加权平均）
    std::vector<Vector> node_sum(n_nodes, Vector(num_strain_components(), 0.0));
    std::vector<Real> node_weight(n_nodes, 0.0);
    
    // 遍历所有单元
    Index gp_offset = 0;
    for (Index elem_id = 0; elem_id < n_elems; ++elem_id) {
        const Element& elem = mesh.element(elem_id);
        ElementType elem_type = elem.type();
        Index n_gp = get_num_gauss_points(elem_type);
        
        // 获取单元节点
        const auto& node_ids = elem.nodes();
        
        // 简化外插：所有高斯点值平均到所有节点
        Vector avg_value(num_strain_components(), 0.0);
        for (Index gp = 0; gp < n_gp; ++gp) {
            const Vector& gp_value = gp_data->get(gp_offset + gp);
            for (Index i = 0; i < num_strain_components(); ++i) {
                avg_value[i] += gp_value[i];
            }
        }
        for (Index i = 0; i < num_strain_components(); ++i) {
            avg_value[i] /= n_gp;
        }
        
        // 分配给所有节点
        for (Index node_id : node_ids) {
            for (Index i = 0; i < num_strain_components(); ++i) {
                node_sum[node_id][i] += avg_value[i];
            }
            node_weight[node_id] += 1.0;
        }
        
        gp_offset += n_gp;
    }
    
    // 计算加权平均
    for (Index node_id = 0; node_id < n_nodes; ++node_id) {
        if (node_weight[node_id] > 0.0) {
            for (Index i = 0; i < num_strain_components(); ++i) {
                node_sum[node_id][i] /= node_weight[node_id];
            }
        }
        node_data->set(node_id, node_sum[node_id]);
    }
    
    FEM_INFO("Extrapolated '" + gp_field + "' to nodes '" + node_field + "'");
}

void PostProcessor::average_to_elements(data::DataManager& data_manager,
                                       const std::string& gp_field,
                                       const std::string& elem_field) {
    // 获取高斯点场
    auto* gp_data = data_manager.get_field<data::VectorData>(gp_field);
    if (!gp_data) {
        FEM_ERROR("PostProcessor: Gauss point field '" + gp_field + "' not found");
        return;
    }
    
    const Mesh& mesh = model_.mesh(0);
    Index n_elems = mesh.num_elements();
    
    // 创建单元场
    auto* elem_data = data_manager.add_field<data::VectorData>(
        elem_field, data::DataLocation::Element, n_elems);
    
    // 遍历所有单元
    Index gp_offset = 0;
    for (Index elem_id = 0; elem_id < n_elems; ++elem_id) {
        ElementType elem_type = mesh.element(elem_id).type();
        Index n_gp = get_num_gauss_points(elem_type);
        
        // 计算平均值
        Vector avg_value(num_strain_components(), 0.0);
        for (Index gp = 0; gp < n_gp; ++gp) {
            const Vector& gp_value = gp_data->get(gp_offset + gp);
            for (Index i = 0; i < num_strain_components(); ++i) {
                avg_value[i] += gp_value[i];
            }
        }
        for (Index i = 0; i < num_strain_components(); ++i) {
            avg_value[i] /= n_gp;
        }
        
        elem_data->set(elem_id, avg_value);
        gp_offset += n_gp;
    }
    
    FEM_INFO("Averaged '" + gp_field + "' to elements '" + elem_field + "'");
}

// ═══════════════════════════════════════════════════════════
// 应力分量提取
// ═══════════════════════════════════════════════════════════

void PostProcessor::extract_stress_component(data::DataManager& data_manager,
                                             const std::string& stress_field,
                                             int component,
                                             const std::string& output_field) {
    // 获取应力场
    auto* stress_data = data_manager.get_field<data::VectorData>(stress_field);
    if (!stress_data) {
        FEM_ERROR("PostProcessor: stress field '" + stress_field + "' not found");
        return;
    }
    
    Index size = stress_data->size();
    data::DataLocation location = stress_data->location();
    
    // 创建标量场
    auto* output_data = data_manager.add_field<data::RealData>(
        output_field, location, size);
    
    // 提取分量
    for (Index i = 0; i < size; ++i) {
        const Vector& stress = stress_data->get(i);
        if (component < stress.size()) {
            output_data->set(i, stress[component]);
        }
    }
    
    FEM_INFO("Extracted component " + std::to_string(component) + 
             " from '" + stress_field + "' to '" + output_field + "'");
}

void PostProcessor::compute_von_mises(data::DataManager& data_manager,
                                     const std::string& stress_field,
                                     const std::string& output_field) {
    // 获取应力场
    auto* stress_data = data_manager.get_field<data::VectorData>(stress_field);
    if (!stress_data) {
        FEM_ERROR("PostProcessor: stress field '" + stress_field + "' not found");
        return;
    }
    
    Index size = stress_data->size();
    data::DataLocation location = stress_data->location();
    
    // 创建标量场
    auto* vm_data = data_manager.add_field<data::RealData>(
        output_field, location, size);
    
    // 计算 von Mises 应力
    for (Index i = 0; i < size; ++i) {
        const Vector& sigma = stress_data->get(i);
        
        Real sigma_xx = sigma[0];
        Real sigma_yy = sigma[1];
        Real sigma_zz = (dimension_ == 3) ? sigma[2] : 0.0;
        Real tau_xy = (dimension_ == 2) ? sigma[2] : sigma[3];
        Real tau_yz = (dimension_ == 3) ? sigma[4] : 0.0;
        Real tau_xz = (dimension_ == 3) ? sigma[5] : 0.0;
        
        // 平均应力
        Real sigma_m = (sigma_xx + sigma_yy + sigma_zz) / 3.0;
        
        // 偏应力
        Real s_xx = sigma_xx - sigma_m;
        Real s_yy = sigma_yy - sigma_m;
        Real s_zz = sigma_zz - sigma_m;
        
        // von Mises: sqrt(1.5 * s:s)
        Real J2 = (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz) / 2.0 +
                  tau_xy * tau_xy + tau_yz * tau_yz + tau_xz * tau_xz;
        
        Real vm = std::sqrt(3.0 * J2);
        vm_data->set(i, vm);
    }
    
    FEM_INFO("Computed von Mises stress '" + output_field + "'");
}

// ═══════════════════════════════════════════════════════════
// 私有辅助函数
// ═══════════════════════════════════════════════════════════

void PostProcessor::compute_element_strain(Index elem_id,
                                          const Mesh& mesh,
                                          const Vector& displacement,
                                          std::vector<Vector>& strains) {
    const Element& elem = mesh.element(elem_id);
    ElementType elem_type = elem.type();
    
    // 获取节点坐标
    const auto& node_ids = elem.nodes();
    std::vector<Vec3> coords;
    for (Index nid : node_ids) {
        coords.push_back(mesh.node(nid).coords());
    }
    
    // 创建形函数对象
    auto shape_func = shape::ShapeFunctionFactory::create(elem_type);
    if (!shape_func) {
        FEM_WARN("PostProcessor: unsupported element type");
        return;
    }
    
    // 获取单元位移向量
    Vector u_elem = get_element_displacement(elem_id, mesh, displacement);
    
    // 获取高斯积分点
    std::vector<Vec3> gauss_points;
    std::vector<Real> weights;
    int order = 2;  // 默认阶数
    shape_func->getGaussPoints(order, gauss_points, weights);
    
    strains.clear();
    strains.reserve(gauss_points.size());
    
    // 在每个高斯点计算应变
    for (const auto& xi : gauss_points) {
        // 计算形函数物理坐标导数
        DenseMatrix dN_dxyz;
        shape_func->computePhysicalDerivatives(xi, coords, dN_dxyz);
        
        // 构造 B 矩阵
        DenseMatrix B;
        if (dimension_ == 2) {
            B = build_B_matrix_2D(dN_dxyz);
        } else {
            B = build_B_matrix_3D(dN_dxyz);
        }
        
        // 应变 = B * u_elem
        Vector strain = B * u_elem;
        strains.push_back(strain);
    }
}

void PostProcessor::compute_stress_at_gp(const Vector& strain,
                                        Vector& stress,
                                        constitutive::StateVariables& state) {
    // 使用材料本构计算应力
    DenseMatrix D;
    material_->computeTangent(strain, D, state);
    
    // σ = D * ε  
    // stress 需要先 resize 到正确的大小
    stress.resize(D.rows());
    for (Index i = 0; i < D.rows(); ++i) {
        stress[i] = 0.0;
        for (Index j = 0; j < D.cols(); ++j) {
            stress[i] += D(i, j) * strain[j];
        }
    }
}

DenseMatrix PostProcessor::build_B_matrix_2D(const DenseMatrix& dN_dxyz) const {
    Index n_nodes = dN_dxyz.cols();
    
    // B 矩阵 (3 × 2*n_nodes)
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
}

DenseMatrix PostProcessor::build_B_matrix_3D(const DenseMatrix& dN_dxyz) const {
    Index n_nodes = dN_dxyz.cols();
    
    // B 矩阵 (6 × 3*n_nodes)
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

Vector PostProcessor::get_element_displacement(Index elem_id,
                                              const Mesh& mesh,
                                              const Vector& displacement) const {
    const auto& node_ids = mesh.element(elem_id).nodes();
    Index n_nodes = node_ids.size();
    
    Vector u_elem(n_nodes * dofs_per_node_);
    
    for (Index i = 0; i < n_nodes; ++i) {
        Index node_id = node_ids[i];
        for (Index d = 0; d < dofs_per_node_; ++d) {
            Index dof = node_id * dofs_per_node_ + d;
            u_elem[i * dofs_per_node_ + d] = displacement[dof];
        }
    }
    
    return u_elem;
}

Index PostProcessor::get_num_gauss_points(ElementType elem_type) const {
    // 默认高斯点数量
    switch (elem_type) {
        case ElementType::Tri3:  return 3;
        case ElementType::Quad4: return 4;
        case ElementType::Tet4:  return 4;
        case ElementType::Brick8: return 8;
        default: return 1;
    }
}

} // namespace postprocess
} // namespace fem
