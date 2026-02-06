#include "physics/heat_conduction.h"
#include <cmath>

namespace fem {

void heat_stiffness(const Vec3* coords, std::size_t n_nodes, std::size_t /*dofs_per_node*/, Real* Ke, void* ctx) {
    const HeatMaterial* material = static_cast<HeatMaterial*>(ctx);
    Real k = material->conductivity;

    // 假设为三角形单元 (3 节点)
    if (n_nodes != 3) {
        // 对于其他单元类型，这里需要扩展
        for (std::size_t i = 0; i < n_nodes * n_nodes; ++i) {
            Ke[i] = 0.0;
        }
        return;
    }

    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    Real area = std::abs(detJ) * 0.5;

    // 参考梯度 dN/dxi:  N0=(-1,-1), N1=(1,0), N2=(0,1)
    // 物理梯度 dN/dx = J^{-T} * dN/dxi
    //   J^{-T} = (1/detJ) * [[ y2-y0, -(y1-y0)],
    //                         [-(x2-x0),  x1-x0]]
    Real dNdxi[3][2] = {{-1, -1}, {1, 0}, {0, 1}};
    Real dNdx[3][2];   // [node][x,y]

    for (int i = 0; i < 3; ++i) {
        dNdx[i][0] = ( (y2-y0)*dNdxi[i][0] - (y1-y0)*dNdxi[i][1] ) / detJ;
        dNdx[i][1] = (-(x2-x0)*dNdxi[i][0] + (x1-x0)*dNdxi[i][1] ) / detJ;
    }

    // Ke[i*3+j] = k * area * dot(dNdx[i], dNdx[j])
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke[i*3+j] = k * area * (dNdx[i][0]*dNdx[j][0] + dNdx[i][1]*dNdx[j][1]);
        }
    }

    // 如果不是3节点，则填充其余部分为0
    for (std::size_t i = 3; i < n_nodes; ++i) {
        for (std::size_t j = 0; j < n_nodes; ++j) {
            Ke[i*n_nodes+j] = 0.0;
        }
        for (std::size_t j = 0; j < 3; ++j) {
            Ke[i*n_nodes+j] = 0.0;
            Ke[j*n_nodes+i] = 0.0;
        }
    }
}

void heat_load(const Vec3* coords, std::size_t n_nodes, std::size_t /*dofs_per_node*/, Real* Fe, void* ctx) {
    const HeatMaterial* material = static_cast<HeatMaterial*>(ctx);
    Real Q = material->source;

    // 假设为三角形单元 (3 节点)
    if (n_nodes != 3) {
        for (std::size_t i = 0; i < n_nodes; ++i) {
            Fe[i] = 0.0;
        }
        return;
    }

    Real x0 = coords[0][0], y0 = coords[0][1];
    Real x1 = coords[1][0], y1 = coords[1][1];
    Real x2 = coords[2][0], y2 = coords[2][1];

    Real detJ = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    Real area = std::abs(detJ) * 0.5;

    for (std::size_t i = 0; i < 3; ++i) {
        Fe[i] = Q * area / 3.0;
    }

    // 如果不是3节点，则填充其余部分为0
    for (std::size_t i = 3; i < n_nodes; ++i) {
        Fe[i] = 0.0;
    }
}

}  // namespace fem