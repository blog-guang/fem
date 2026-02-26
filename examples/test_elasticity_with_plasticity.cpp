/**
 * test_elasticity_with_plasticity.cpp - 演示如何使用自定义材料模型
 * 
 * 展示 ElasticityUnified 集成 J2 塑性材料
 */

#include "physics/elasticity_unified.h"
#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "solver/cg.h"
#include "core/logger.h"
#include <iostream>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "=== ElasticityUnified with Custom Materials ===\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 示例 1: 2D 平面应力弹性材料
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 1: 2D Plane Stress (IsotropicElastic)\n";
    {
        Model model("plane_stress");
        int mat_id = model.add_material("steel");
        int mesh_id = model.add_mesh("mesh", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        MeshGenerator::generate_unit_square_tri(2, 2, mesh);
        
        Real E = 200e9;   // 200 GPa
        Real nu = 0.3;
        
        // 创建 2D 平面应力材料
        IsotropicElastic material(E, nu, 2, true);  // 2D, plane_stress
        ElasticityUnified physics(&material, 2);
        
        std::cout << "  Material: IsotropicElastic\n";
        std::cout << "  E  = " << E / 1e9 << " GPa\n";
        std::cout << "  nu = " << nu << "\n";
        std::cout << "  mu = " << material.mu() / 1e9 << " GPa\n";
        std::cout << "  Dimension: " << physics.dimension() << "D\n";
        
        // 装配
        Assembler assembler(model, 2);
        auto elem_func = [&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            physics.compute_element(elem_id, m, Ke, Fe);
        };
        assembler.assemble(elem_func);
        
        std::cout << "  Global DOFs: " << assembler.num_dofs() << "\n";
        std::cout << "  ✓ Assembly successful\n\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 2: 高级接口（使用 IsotropicElastic 对象）
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 2: Advanced interface (IsotropicElastic object)\n";
    {
        Model model("advanced_elastic");
        int mat_id = model.add_material("aluminum");
        int mesh_id = model.add_mesh("mesh", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        MeshGenerator::generate_unit_square_quad(2, 2, mesh);
        
        Real E = 70e9;    // 70 GPa (铝)
        Real nu = 0.33;
        
        // 创建材料对象
        IsotropicElastic material(E, nu, 2, true);  // 2D, plane_stress
        
        // 高级构造：传入材料指针
        ElasticityUnified physics(&material, 2);
        
        std::cout << "  Material: IsotropicElastic\n";
        std::cout << "  E  = " << E / 1e9 << " GPa\n";
        std::cout << "  mu = " << material.mu() / 1e9 << " GPa\n";
        std::cout << "  Dimension: " << physics.dimension() << "D\n";
        
        // 装配
        Assembler assembler(model, 2);
        auto elem_func = [&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            physics.compute_element(elem_id, m, Ke, Fe);
        };
        assembler.assemble(elem_func);
        
        std::cout << "  Global DOFs: " << assembler.num_dofs() << "\n";
        std::cout << "  ✓ Assembly successful\n\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 3: J2 塑性材料（弹塑性）
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 3: J2 Plasticity material (elasto-plastic)\n";
    {
        Model model("plasticity");
        int mat_id = model.add_material("mild_steel");
        int mesh_id = model.add_mesh("mesh", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        MeshGenerator::generate_unit_square_tri(3, 3, mesh);
        
        Real E = 200e9;       // 200 GPa
        Real nu = 0.3;
        Real sigma_y = 250e6; // 250 MPa (屈服应力)
        Real H = 2e9;         // 2 GPa (硬化模量)
        
        // 创建 J2 塑性材料
        J2Plasticity material(E, nu, sigma_y, H, 2);  // 2D
        
        // 使用塑性材料
        ElasticityUnified physics(&material, 2);
        
        std::cout << "  Material: J2Plasticity\n";
        std::cout << "  E       = " << E / 1e9 << " GPa\n";
        std::cout << "  sigma_y = " << sigma_y / 1e6 << " MPa\n";
        std::cout << "  H       = " << H / 1e9 << " GPa\n";
        std::cout << "  Dimension: " << physics.dimension() << "D\n";
        
        // 装配（对于塑性材料，这只是初始弹性刚度）
        Assembler assembler(model, 2);
        auto elem_func = [&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            physics.compute_element(elem_id, m, Ke, Fe);
        };
        assembler.assemble(elem_func);
        
        std::cout << "  Global DOFs: " << assembler.num_dofs() << "\n";
        std::cout << "  ✓ Initial elastic assembly successful\n";
        std::cout << "  Note: For plasticity, use Newton-Raphson with tangent updates\n\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 4: 3D 弹塑性
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 4: 3D J2 Plasticity\n";
    {
        Model model("plasticity_3d");
        int mat_id = model.add_material("steel_3d");
        int mesh_id = model.add_mesh("mesh", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        MeshGenerator::generate_unit_cube_tet(2, 2, 2, mesh);
        
        Real E = 200e9;
        Real nu = 0.3;
        Real sigma_y = 250e6;
        Real H = 2e9;
        
        // 3D J2 塑性
        J2Plasticity material(E, nu, sigma_y, H, 3);  // 3D
        ElasticityUnified physics(&material, 3);
        
        std::cout << "  Material: J2Plasticity (3D)\n";
        std::cout << "  Dimension: " << physics.dimension() << "D\n";
        std::cout << "  Nodes: " << mesh.num_nodes() << "\n";
        std::cout << "  Elements: " << mesh.num_elements() << "\n";
        
        Assembler assembler(model, 3);
        auto elem_func = [&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            physics.compute_element(elem_id, m, Ke, Fe);
        };
        assembler.assemble(elem_func);
        
        std::cout << "  Global DOFs: " << assembler.num_dofs() << "\n";
        std::cout << "  ✓ Assembly successful\n\n";
    }
    
    std::cout << "=== All examples completed successfully! ===\n";
    
    return 0;
}
