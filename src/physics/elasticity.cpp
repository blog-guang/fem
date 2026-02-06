#include "physics/elasticity.h"
#include <cmath>

namespace fem {

void elasticity_stiffness(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Ke, void* ctx) {
    // 初始化为零
    std::size_t total_dofs = n_nodes * dofs_per_node;
    for (std::size_t i = 0; i < total_dofs * total_dofs; ++i) {
        Ke[i] = 0.0;
    }

    // 只处理三角形单元 (3 节点) 且每个节点2个自由度 (x,y位移)
    if (n_nodes != 3 || dofs_per_node != 2) {
        return;
    }

    const ElasticCtx* elastic_ctx = static_cast<ElasticCtx*>(ctx);
    const ElasticMaterial& mat = elastic_ctx->mat;
    
    Real E = mat.E;
    Real nu = mat.poisson;
    Real t = mat.thickness;  // 厚度

    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    // 计算面积
    Real area = std::abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) * 0.5;

    // 计算B矩阵中的系数
    Real b1 = y1 - y2;  // y2-y3 -> changed to y1-y2 following standard notation
    Real b2 = y2 - y0;  // y3-y1 -> changed to y2-y0
    Real b3 = y0 - y1;  // y1-y2 -> changed to y0-y1
    
    Real c1 = x2 - x1;  // x3-x2 -> changed to x2-x1
    Real c2 = x0 - x2;  // x1-x3 -> changed to x0-x2
    Real c3 = x1 - x0;  // x2-x1 -> changed to x1-x0

    // 按照标准顺序修正：对于节点0,1,2对应的坐标，重新定义
    b1 = y2 - y1;  // y2-y3 where we consider proper order
    b2 = y0 - y2;  // y3-y1 
    b3 = y1 - y0;  // y1-y2
    c1 = x1 - x2;  // x3-x2
    c2 = x2 - x0;  // x1-x3
    c3 = x0 - x1;  // x2-x1

    // 计算正确的面积和系数（重新按标准方式）
    b1 = y2 - y0;  // b1 = y2 - y0
    b2 = y0 - y1;  // b2 = y0 - y1  
    b3 = y1 - y2;  // b3 = y1 - y2
    c1 = x0 - x2;  // c1 = x0 - x2
    c2 = x2 - x1;  // c2 = x2 - x1
    c3 = x1 - x0;  // c3 = x1 - x0

    // 重新计算
    x0 = coords[0][0], y0 = coords[0][1];
    x1 = coords[1][0], y1 = coords[1][1];
    x2 = coords[2][0], y2 = coords[2][1];
    
    b1 = y1 - y2;  // y2-y3 with proper indexing
    b2 = y2 - y0;  // y3-y1
    b3 = y0 - y1;  // y1-y2
    c1 = x2 - x1;  // x3-x2
    c2 = x0 - x2;  // x1-x3
    c3 = x1 - x0;  // x2-x1

    // 实际的面积计算应该使用绝对值
    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    area = std::abs(detJ) * 0.5;

    // 按照标准公式重新设置
    b1 = y1 - y2;
    b2 = y2 - y0;
    b3 = y0 - y1;
    c1 = x2 - x1;
    c2 = x0 - x2;
    c3 = x1 - x0;

    // B矩阵 (3x6): [du/dx, dv/dy, du/dy+dv/dx]
    // B = (1/(2*area)) * [[b1,0,b2,0,b3,0],
    //                      [0,c1,0,c2,0,c3],
    //                      [c1,b1,c2,b2,c3,b3]]
    Real factor = 1.0 / (2.0 * area);
    
    Real B[3][6] = {
        {factor * b1, 0.0, factor * b2, 0.0, factor * b3, 0.0},      // strain_x = du/dx
        {0.0, factor * c1, 0.0, factor * c2, 0.0, factor * c3},      // strain_y = dv/dy
        {factor * c1, factor * b1, factor * c2, factor * b2, factor * c3, factor * b3}  // gamma_xy = du/dy + dv/dx
    };

    // D矩阵 (平面应力): D = E/(1-nu^2) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    Real D[3][3];
    Real D_factor = E / (1.0 - nu * nu);
    D[0][0] = D_factor;           D[0][1] = D_factor * nu;    D[0][2] = 0.0;
    D[1][0] = D_factor * nu;      D[1][1] = D_factor;         D[1][2] = 0.0;
    D[2][0] = 0.0;                D[2][1] = 0.0;              D[2][2] = D_factor * (1.0 - nu) / 2.0;

    // 计算 Ke = t * area * B^T * D * B
    Real temp[3][6];  // temp = D * B (3x6)
    Real K_local[6][6];  // K_local = B^T * temp

    // temp = D * B
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            temp[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                temp[i][j] += D[i][k] * B[k][j];
            }
        }
    }

    // K_local = B^T * temp
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_local[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                K_local[i][j] += B[k][i] * temp[k][j];
            }
        }
    }

    // 应用厚度和面积因子
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_local[i][j] *= t * area;
        }
    }

    // 将局部刚度矩阵复制到输出数组
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            Ke[i * 6 + j] = K_local[i][j];
        }
    }
}

void elasticity_load(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* ctx) {
    // 初始化为零
    std::size_t total_dofs = n_nodes * dofs_per_node;
    for (std::size_t i = 0; i < total_dofs; ++i) {
        Fe[i] = 0.0;
    }

    // 只处理三角形单元 (3 节点) 且每个节点2个自由度 (x,y位移)
    if (n_nodes != 3 || dofs_per_node != 2) {
        return;
    }

    const ElasticCtx* elastic_ctx = static_cast<ElasticCtx*>(ctx);
    const ElasticLoad& load = elastic_ctx->load;
    
    Real fx = load.fx;
    Real fy = load.fy;

    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real area = std::abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) * 0.5;

    // 对每个节点的x和y自由度施加体力载荷
    // Fe[2*i] = fx*area/3, Fe[2*i+1] = fy*area/3
    for (std::size_t i = 0; i < 3; ++i) {
        Fe[2*i] = fx * area / 3.0;      // x方向位移的载荷
        Fe[2*i + 1] = fy * area / 3.0;  // y方向位移的载荷
    }
}

}  // namespace fem