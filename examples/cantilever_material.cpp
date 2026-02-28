/**
 * cantilever_material.cpp - 使用材料本构模型的悬臂梁分析
 * 
 * 演示：
 * 1. IsotropicElastic - 小载荷（线性）
 * 2. J2Plasticity - 大载荷（非线性，塑性）
 * 3. 增量加载法
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_material.h"
#include "material/isotropic_elastic.h"
#include "material/j2_plasticity.h"
#include "math/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include "core/logger.h"
#include <memory>
#include <iomanip>
#include <iostream>

using namespace fem;
using namespace fem::constitutive;
using namespace fem::physics;

void cantilever_elastic(Real load_factor = -10.0) {
    FEM_INFO("=== 悬臂梁：IsotropicElastic (线性分析) ===");
    
    // 创建模型
    Model model("CantileverElastic");
    
    int mat_id = model.add_material("Steel");
    model.material(mat_id).set_property("E", 200e3);
    model.material(mat_id).set_property("nu", 0.3);
    
    int mesh_id = model.add_mesh("beam", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成网格：20x5 (长4 x 高1)
    MeshGenerator::generate_unit_square_tri(20, 5, mesh);
    
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        auto& coords = mesh.node(i).coords();
        coords[0] *= 4.0;  // 长度4
        coords[1] *= 1.0;  // 高度1
    }
    
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("网格: " + std::to_string(mesh.num_nodes()) + " 节点, " +
             std::to_string(mesh.num_elements()) + " 单元");
    
    // 创建材料本构
    Real E = 200e3, nu = 0.3;
    auto material = std::make_shared<IsotropicElastic>(E, nu, 2, true);
    FEM_INFO("材料: " + material->typeName());
    
    // 创建物理模块
    ElasticityWithMaterial elast_mat(material, 0.1);  // 厚度0.1
    elast_mat.initialize(mesh.num_elements());
    
    // 装配
    Timer timer;
    timer.start();
    
    Assembler assembler(model, 2);
    
    auto elem_func = [&](Index elem_id, const Mesh& m,
                        DenseMatrix& Ke, Vector& Fe) {
        elast_mat.compute_stiffness(elem_id, m, Ke);
        Fe.resize(Ke.rows(), 0.0);
    };
    
    assembler.assemble(elem_func);
    
    FEM_INFO("装配耗时: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 边界条件：左边固定
    std::vector<DirichletBC> dirichlet_bcs = {
        {"left", 0, 0.0},
        {"left", 1, 0.0}
    };
    assembler.apply_dirichlet(dirichlet_bcs);
    
    // 载荷：顶部均布
    std::vector<NeumannBC> neumann_bcs = {
        {"top", 1, load_factor}
    };
    assembler.apply_neumann(neumann_bcs);
    
    // 求解
    timer.start();
    
    CGSolver solver;
    solver.set_tolerance(1e-10);
    
    Vector u(assembler.num_dofs(), 0.0);
    auto result = solver.solve(assembler.matrix(), assembler.rhs(), u);
    
    FEM_INFO("求解: " + std::to_string(result.iterations) + " 次迭代, " +
             std::to_string(timer.elapsed_s()) + "s");
    
    // 更新应力
    for (Index elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        Vector Fe(6);
        elast_mat.update_stress(elem_id, mesh, u, Fe);
    }
    
    // 后处理：最大位移
    Real max_uy = 0.0;
    Index max_node = 0;
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Real uy = std::abs(u[2*i + 1]);
        if (uy > max_uy) {
            max_uy = uy;
            max_node = i;
        }
    }
    
    const auto& coords = mesh.node(max_node).coords();
    FEM_INFO("最大位移 u_y = " + std::to_string(max_uy) +
             " at (" + std::to_string(coords[0]) + ", " + std::to_string(coords[1]) + ")");
    
    // 输出VTK
    VTKWriter vtk("cantilever_elastic.vtu");
    vtk.write_mesh(mesh);
    vtk.write_point_data("displacement", u, 2);
    
    FEM_INFO("输出: cantilever_elastic.vtu");
}

void cantilever_plastic(Real load_factor = -100.0, int num_steps = 5) {
    FEM_INFO("\n=== 悬臂梁：J2Plasticity (非线性塑性分析) ===");
    
    Model model("CantileverPlastic");
    
    int mat_id = model.add_material("MildSteel");
    model.material(mat_id).set_property("E", 200e3);
    model.material(mat_id).set_property("nu", 0.3);
    
    int mesh_id = model.add_mesh("beam", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_tri(20, 5, mesh);
    
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        auto& coords = mesh.node(i).coords();
        coords[0] *= 4.0;
        coords[1] *= 1.0;
    }
    
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("网格: " + std::to_string(mesh.num_nodes()) + " 节点, " +
             std::to_string(mesh.num_elements()) + " 单元");
    
    // 创建塑性材料
    Real E = 200e3, nu = 0.3, sigma_y = 250.0, H = 5000.0;
    auto material = std::make_shared<J2Plasticity>(E, nu, sigma_y, H, 2);
    FEM_INFO("材料: " + material->typeName());
    FEM_INFO("  σ_y0 = " + std::to_string(sigma_y) + " MPa");
    FEM_INFO("  H = " + std::to_string(H) + " MPa");
    
    ElasticityWithMaterial elast_mat(material, 0.1);
    elast_mat.initialize(mesh.num_elements());
    
    // 增量加载
    Vector u_total(2 * mesh.num_nodes(), 0.0);
    
    FEM_INFO("\n增量加载: " + std::to_string(num_steps) + " 步");
    
    Timer total_timer;
    total_timer.start();
    
    for (int step = 1; step <= num_steps; ++step) {
        FEM_INFO("\n--- 加载步 " + std::to_string(step) + "/" + std::to_string(num_steps) + " ---");
        
        Timer step_timer;
        step_timer.start();
        
        // 装配
        Assembler assembler(model, 2);
        
        auto elem_func = [&](Index elem_id, const Mesh& m,
                            DenseMatrix& Ke, Vector& Fe) {
            elast_mat.compute_stiffness(elem_id, m, Ke);
            Fe.resize(Ke.rows(), 0.0);
        };
        
        assembler.assemble(elem_func);
        
        // 边界条件
        std::vector<DirichletBC> dirichlet_bcs = {
            {"left", 0, 0.0},
            {"left", 1, 0.0}
        };
        assembler.apply_dirichlet(dirichlet_bcs);
        
        // 增量载荷
        std::vector<NeumannBC> neumann_bcs = {
            {"top", 1, load_factor / num_steps}
        };
        assembler.apply_neumann(neumann_bcs);
        
        // 求解位移增量
        CGSolver solver;
        solver.set_tolerance(1e-10);
        
        Vector u_inc(assembler.num_dofs(), 0.0);
        auto result = solver.solve(assembler.matrix(), assembler.rhs(), u_inc);
        
        FEM_INFO("  求解: " + std::to_string(result.iterations) + " 次迭代");
        
        // 更新应力和状态
        for (Index elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
            Vector Fe(6);
            elast_mat.update_stress(elem_id, mesh, u_inc, Fe);
        }
        
        // 累积位移
        u_total += u_inc;
        
        FEM_INFO("  步骤耗时: " + std::to_string(step_timer.elapsed_s()) + "s");
    }
    
    FEM_INFO("\n总耗时: " + std::to_string(total_timer.elapsed_s()) + "s");
    
    // 统计塑性状态
    int plastic_elements = 0;
    Real max_plastic_strain = 0.0;
    
    for (Index elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        const auto& state = elast_mat.get_state(elem_id);
        if (state.equiv_plastic_strain > 1e-10) {
            plastic_elements++;
            max_plastic_strain = std::max(max_plastic_strain, 
                                         state.equiv_plastic_strain);
        }
    }
    
    FEM_INFO("\n塑性统计:");
    FEM_INFO("  塑性单元: " + std::to_string(plastic_elements) + " / " +
             std::to_string(mesh.num_elements()) +
             " (" + std::to_string(100.0 * plastic_elements / mesh.num_elements()) + "%)");
    FEM_INFO("  最大等效塑性应变: " + std::to_string(max_plastic_strain));
    
    // 最大位移
    Real max_uy = 0.0;
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        max_uy = std::max(max_uy, std::abs(u_total[2*i + 1]));
    }
    FEM_INFO("  最大位移 u_y = " + std::to_string(max_uy));
    
    // 输出VTK
    VTKWriter vtk("cantilever_plastic.vtu");
    vtk.write_mesh(mesh);
    vtk.write_point_data("displacement", u_total, 2);
    
    // 输出塑性应变
    Vector plastic_strain(mesh.num_elements());
    for (Index elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        plastic_strain[elem_id] = elast_mat.get_state(elem_id).equiv_plastic_strain;
    }
    vtk.write_cell_data("plastic_strain", plastic_strain);
    
    FEM_INFO("输出: cantilever_plastic.vtu");
}

int main() {
    std::cout << "\n";
    std::cout << "███████████████████████████████████████████████████████████\n";
    std::cout << "  悬臂梁材料本构分析\n";
    std::cout << "  Author: Math Agent 🧮\n";
    std::cout << "███████████████████████████████████████████████████████████\n\n";
    
    try {
        // 测试1：弹性分析
        cantilever_elastic(-10.0);
        
        // 测试2：塑性分析（大载荷）
        cantilever_plastic(-100.0, 5);
        
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "  分析完成！\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
