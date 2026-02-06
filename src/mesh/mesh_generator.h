#pragma once

#include "mesh/mesh.h"

namespace fem {

// ── 结构化网格生成器 ──
// 生成规则网格，用于快速原型和测试

// 单位正方形 [0,1]x[0,1], 三角形网格
// nx, ny: 每个方向的划分数
// 每个方格分成 2 个三角形
[[nodiscard]] Mesh generate_unit_square_tri(std::size_t nx, std::size_t ny);

// 单位正方形 [0,1]x[0,1], 四边形网格
[[nodiscard]] Mesh generate_unit_square_quad(std::size_t nx, std::size_t ny);

}  // namespace fem
