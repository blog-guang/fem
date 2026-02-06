#pragma once

#include "mesh/mesh.h"
#include "core/types.h"

namespace fem {

/**
 * 网格生成器
 * 
 * 生成基于新 Mesh 架构的结构化网格。
 */
class MeshGenerator {
public:
    /**
     * 生成 2D 单位正方形三角形网格 [0,1]x[0,1]
     * 
     * @param nx x 方向单元数
     * @param ny y 方向单元数
     * @param mesh 目标 Mesh (已绑定材料)
     */
    static void generate_unit_square_tri(int nx, int ny, Mesh& mesh);
    
    /**
     * 生成 2D 单位正方形四边形网格 [0,1]x[0,1]
     * 
     * @param nx x 方向单元数
     * @param ny y 方向单元数
     * @param mesh 目标 Mesh (已绑定材料)
     */
    static void generate_unit_square_quad(int nx, int ny, Mesh& mesh);
    
    /**
     * 生成 3D 单位立方体四面体网格 [0,1]x[0,1]x[0,1]
     * 
     * @param nx x 方向单元数
     * @param ny y 方向单元数
     * @param nz z 方向单元数
     * @param mesh 目标 Mesh (已绑定材料)
     */
    static void generate_unit_cube_tet(int nx, int ny, int nz, Mesh& mesh);
    
    /**
     * 生成 3D 单位立方体六面体网格 [0,1]x[0,1]x[0,1]
     * 
     * @param nx x 方向单元数
     * @param ny y 方向单元数
     * @param nz z 方向单元数
     * @param mesh 目标 Mesh (已绑定材料)
     */
    static void generate_unit_cube_brick(int nx, int ny, int nz, Mesh& mesh);
    
    /**
     * 自动识别边界并命名
     * 
     * 2D: "left", "right", "bottom", "top"
     * 3D: 加上 "front", "back"
     */
    static void identify_boundaries_2d(Mesh& mesh);
    static void identify_boundaries_3d(Mesh& mesh);
};

}  // namespace fem
