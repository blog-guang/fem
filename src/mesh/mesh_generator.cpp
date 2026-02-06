#include "mesh/mesh_generator.h"
#include "core/logger.h"

namespace fem {

//
// 节点编号示意 (nx=3, ny=2):
//
//   6 ── 7 ── 8 ── 9
//   |    |    |    |
//   3 ── 4 ── 5 ── 6   ← 每行 (nx+1) 个节点
//   |    |    |    |
//   0 ── 1 ── 2 ── 3
//
// node(ix, iy) = iy * (nx+1) + ix
//

Mesh generate_unit_square_tri(std::size_t nx, std::size_t ny) {
    Mesh mesh;

    Real dx = 1.0 / static_cast<Real>(nx);
    Real dy = 1.0 / static_cast<Real>(ny);

    // ── 添加节点 ──
    for (std::size_t iy = 0; iy <= ny; ++iy) {
        for (std::size_t ix = 0; ix <= nx; ++ix) {
            mesh.add_node({ix * dx, iy * dy, 0.0});
        }
    }

    // ── 添加单元 (每个方格 → 2 个三角形) ──
    //   (ix,iy+1) ── (ix+1,iy+1)
    //      |  \  2  |
    //      | 1  \   |
    //   (ix,iy) ── (ix+1,iy)
    //
    auto nid = [&](std::size_t ix, std::size_t iy) -> Index {
        return iy * (nx + 1) + ix;
    };

    for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
            Index n00 = nid(ix,   iy);
            Index n10 = nid(ix+1, iy);
            Index n01 = nid(ix,   iy+1);
            Index n11 = nid(ix+1, iy+1);

            // 三角形 1: 下左
            Index tri1[] = {n00, n10, n01};
            mesh.add_cell(ElementType::Triangle2D, tri1, 3);

            // 三角形 2: 上右
            Index tri2[] = {n10, n11, n01};
            mesh.add_cell(ElementType::Triangle2D, tri2, 3);
        }
    }

    // ── 边界标记 ──
    std::vector<Index> bottom, top, left, right;
    for (std::size_t ix = 0; ix <= nx; ++ix) {
        bottom.push_back(nid(ix, 0));
        top.push_back(nid(ix, ny));
    }
    for (std::size_t iy = 0; iy <= ny; ++iy) {
        left.push_back(nid(0,  iy));
        right.push_back(nid(nx, iy));
    }
    mesh.add_boundary("bottom", bottom);
    mesh.add_boundary("top",    top);
    mesh.add_boundary("left",   left);
    mesh.add_boundary("right",  right);

    FEM_INFO("Generated unit_square_tri: " +
             std::to_string(nx) + "x" + std::to_string(ny));
    mesh.print_info();
    return mesh;
}

Mesh generate_unit_square_quad(std::size_t nx, std::size_t ny) {
    Mesh mesh;

    Real dx = 1.0 / static_cast<Real>(nx);
    Real dy = 1.0 / static_cast<Real>(ny);

    // ── 添加节点 (同上) ──
    for (std::size_t iy = 0; iy <= ny; ++iy) {
        for (std::size_t ix = 0; ix <= nx; ++ix) {
            mesh.add_node({ix * dx, iy * dy, 0.0});
        }
    }

    auto nid = [&](std::size_t ix, std::size_t iy) -> Index {
        return iy * (nx + 1) + ix;
    };

    // ── 添加单元 (每个方格 → 1 个四边形, 逆时针) ──
    for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
            Index quad[] = {nid(ix, iy), nid(ix+1, iy), nid(ix+1, iy+1), nid(ix, iy+1)};
            mesh.add_cell(ElementType::Quad2D, quad, 4);
        }
    }

    // ── 边界标记 (同上) ──
    std::vector<Index> bottom, top, left, right;
    for (std::size_t ix = 0; ix <= nx; ++ix) {
        bottom.push_back(nid(ix, 0));
        top.push_back(nid(ix, ny));
    }
    for (std::size_t iy = 0; iy <= ny; ++iy) {
        left.push_back(nid(0,  iy));
        right.push_back(nid(nx, iy));
    }
    mesh.add_boundary("bottom", bottom);
    mesh.add_boundary("top",    top);
    mesh.add_boundary("left",   left);
    mesh.add_boundary("right",  right);

    FEM_INFO("Generated unit_square_quad: " +
             std::to_string(nx) + "x" + std::to_string(ny));
    mesh.print_info();
    return mesh;
}

}  // namespace fem
