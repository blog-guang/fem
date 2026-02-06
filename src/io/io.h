#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include <string>
#include <vector>

namespace fem {

// ── VTK Legacy 写入器 ──
class VTKWriter {
public:
    explicit VTKWriter(const std::string& output_dir);

    // 写入 .vtk 文件
    // fields: 节点标量场列表 {name, values}
    struct ScalarField {
        std::string name;
        const Vector* values;  // size = num_nodes
    };

    void write(const std::string& filename,
               const Mesh&        mesh,
               const std::vector<ScalarField>& fields = {}) const;

private:
    std::string output_dir_;
};

}  // namespace fem
