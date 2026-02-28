/**
 * test_data_manager.cpp - 数据管理系统示例
 * 
 * 展示如何使用 DataManager 和 FieldData 管理网格数据
 */

#include "data/data_manager.h"
#include "data/field_data.h"
#include "core/logger.h"
#include <iostream>

using namespace fem;
using namespace fem::data;

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "=== FEM Data Management System Demo ===\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 示例 1: 基本的场数据使用
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 1: Basic Field Data Usage\n";
    std::cout << "──────────────────────────────────\n";
    {
        Index n_nodes = 100;
        
        // 创建标量场（温度）
        RealData temperature("temperature", DataLocation::Node, n_nodes, 20.0);
        
        // 设置部分节点温度
        for (Index i = 0; i < 10; ++i) {
            temperature.set(i, 100.0 + i * 5.0);
        }
        
        // 读取温度
        std::cout << "Temperature at node 0: " << temperature.get(0) << " K\n";
        std::cout << "Temperature at node 5: " << temperature[5] << " K\n";
        std::cout << "Temperature at node 50: " << temperature[50] << " K (default)\n";
        
        temperature.print_info();
        std::cout << "\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 2: 矢量场和张量场
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 2: Vector and Tensor Fields\n";
    std::cout << "────────────────────────────────────\n";
    {
        Index n_nodes = 50;
        
        // 创建矢量场（位移）
        VectorData displacement("displacement", DataLocation::Node, n_nodes);
        
        // 设置位移
        Vector u(3);
        u[0] = 0.01;  // u_x
        u[1] = 0.02;  // u_y
        u[2] = 0.00;  // u_z
        displacement.set(10, u);
        
        // 读取位移
        const Vector& u_read = displacement.get(10);
        std::cout << "Displacement at node 10: [" 
                  << u_read[0] << ", " 
                  << u_read[1] << ", " 
                  << u_read[2] << "]\n";
        
        // 创建张量场（应力）
        Index n_elems = 40;
        Index n_gp = 4;  // 每单元 4 个高斯点
        TensorData stress("stress", DataLocation::GaussPoint, n_elems * n_gp);
        
        // 设置应力张量（使用 Voigt 记号：6 个分量）
        DenseMatrix sigma(6, 1);  // [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]
        sigma(0, 0) = 100.0;  // σ_xx
        sigma(1, 0) = 50.0;   // σ_yy
        sigma(2, 0) = 30.0;   // σ_zz
        sigma(3, 0) = 10.0;   // τ_xy
        sigma(4, 0) = 0.0;    // τ_yz
        sigma(5, 0) = 0.0;    // τ_xz
        
        Index elem_id = 5;
        Index gp_idx = 2;
        stress.set_gauss_point(elem_id, gp_idx, n_gp, sigma);
        
        // 读取应力
        const DenseMatrix& sigma_read = stress.get_gauss_point(elem_id, gp_idx, n_gp);
        std::cout << "Stress at elem 5, GP 2: σ_xx = " << sigma_read(0, 0) << " MPa\n";
        
        std::cout << "\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 3: DataManager 使用
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 3: DataManager Usage\n";
    std::cout << "─────────────────────────────\n";
    {
        DataManager manager;
        
        Index n_nodes = 100;
        Index n_elems = 80;
        Index n_gp_per_elem = 4;
        
        // 注册节点场
        manager.add_field<RealData>("temperature", DataLocation::Node, n_nodes, 293.15);
        manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes);
        manager.add_field<VectorData>("velocity", DataLocation::Node, n_nodes);
        
        // 注册单元场
        manager.add_field<IntData>("material_id", DataLocation::Element, n_elems, 1);
        manager.add_field<RealData>("element_volume", DataLocation::Element, n_elems);
        
        // 注册高斯点场
        manager.add_field<TensorData>("stress", DataLocation::GaussPoint, n_elems * n_gp_per_elem);
        manager.add_field<TensorData>("strain", DataLocation::GaussPoint, n_elems * n_gp_per_elem);
        
        // 访问场数据
        auto* temp = manager.get_field<RealData>("temperature");
        temp->set(10, 350.0);
        
        auto* disp = manager.get_field<VectorData>("displacement");
        Vector u(3, 0.0);
        u[0] = 0.01;
        disp->set(10, u);
        
        auto* mat_id = manager.get_field<IntData>("material_id");
        mat_id->set(5, 2);  // 单元 5 使用材料 2
        
        // 打印统计信息
        std::cout << "Total fields: " << manager.num_fields() << "\n";
        std::cout << "Node fields: " << manager.num_fields(DataLocation::Node) << "\n";
        std::cout << "Element fields: " << manager.num_fields(DataLocation::Element) << "\n";
        std::cout << "GaussPoint fields: " << manager.num_fields(DataLocation::GaussPoint) << "\n";
        
        std::cout << "\nField names:\n";
        for (const auto& name : manager.get_field_names()) {
            std::cout << "  - " << name << "\n";
        }
        
        std::cout << "\n";
        manager.print_info();
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 4: 网格尺寸变化处理
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "\nExample 4: Mesh Refinement Handling\n";
    std::cout << "────────────────────────────────────\n";
    {
        DataManager manager;
        
        Index n_nodes_old = 100;
        
        // 初始网格数据
        manager.add_field<RealData>("temperature", DataLocation::Node, n_nodes_old, 20.0);
        manager.add_field<VectorData>("displacement", DataLocation::Node, n_nodes_old);
        
        std::cout << "Initial mesh: " << n_nodes_old << " nodes\n";
        
        // 网格细化后
        Index n_nodes_new = 250;
        std::cout << "After refinement: " << n_nodes_new << " nodes\n";
        
        // 批量调整所有节点场的大小
        manager.resize_all(DataLocation::Node, n_nodes_new);
        
        auto* temp = manager.get_field<RealData>("temperature");
        std::cout << "Temperature field size: " << temp->size() << "\n";
        
        std::cout << "\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 示例 5: 整型和布尔场
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Example 5: Integer and Boolean Fields\n";
    std::cout << "──────────────────────────────────────\n";
    {
        Index n_nodes = 50;
        Index n_elems = 40;
        
        // 整型场：边界条件标记
        IntData bc_marker("bc_marker", DataLocation::Node, n_nodes, 0);
        bc_marker.set(0, 1);   // 节点 0: 固定边界
        bc_marker.set(10, 2);  // 节点 10: 对称边界
        
        std::cout << "BC marker at node 0: " << bc_marker[0] << " (fixed)\n";
        std::cout << "BC marker at node 10: " << bc_marker[10] << " (symmetric)\n";
        std::cout << "BC marker at node 20: " << bc_marker[20] << " (free)\n";
        
        // 布尔场：单元激活状态
        BoolData is_active("is_active", DataLocation::Element, n_elems, true);
        is_active.set(5, false);  // 单元 5 失效
        
        std::cout << "Element 3 active: " << (is_active[3] ? "yes" : "no") << "\n";
        std::cout << "Element 5 active: " << (is_active[5] ? "yes" : "no") << "\n";
        
        std::cout << "\n";
    }
    
    std::cout << "=== All examples completed successfully! ===\n";
    
    return 0;
}
