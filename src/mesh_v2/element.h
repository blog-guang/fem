#pragma once

#include "core/types.h"
#include <vector>
#include <array>
#include <memory>

namespace fem {
namespace v2 {

// ElementType 定义在 core/types.h 中

// ═══ Element 抽象基类 ═══
class Element {
public:
    virtual ~Element() = default;
    
    // 基本属性
    virtual ElementType type() const = 0;
    virtual int dimension() const = 0;      // 0=Node, 1=Edge, 2=Face, 3=Volume
    virtual int num_nodes() const = 0;
    virtual int order() const = 0;          // 1=线性, 2=二次
    
    // 节点访问
    virtual const std::vector<Index>& nodes() const = 0;
    
    // 局部编号 (在 Mesh 内的索引)
    void set_local_id(Index id) { local_id_ = id; }
    Index local_id() const { return local_id_; }
    
    // 标签/属性
    void set_tag(int tag) { tag_ = tag; }
    int tag() const { return tag_; }
    
protected:
    Index local_id_{0};
    int tag_{0};
};

// ═══ Node (0D) ═══
class Node : public Element {
public:
    Node(const Vec3& coords) : coords_(coords), dummy_nodes_() {}
    
    ElementType type() const override { return ElementType::Node; }
    int dimension() const override { return 0; }
    int num_nodes() const override { return 1; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return dummy_nodes_; }
    
    const Vec3& coords() const { return coords_; }
    Vec3& coords() { return coords_; }
    
private:
    Vec3 coords_;
    std::vector<Index> dummy_nodes_;  // Node 本身不引用其他节点
};

// ═══ Edge (1D) ═══
class Edge : public Element {
public:
    int dimension() const override { return 1; }
};

// 2节点线性边
class Edge2 : public Edge {
public:
    Edge2(Index n0, Index n1) : nodes_{n0, n1} {}
    
    ElementType type() const override { return ElementType::Edge2; }
    int num_nodes() const override { return 2; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;
};

// 3节点二次边
class Edge3 : public Edge {
public:
    Edge3(Index n0, Index n1, Index n_mid) : nodes_{n0, n1, n_mid} {}
    
    ElementType type() const override { return ElementType::Edge3; }
    int num_nodes() const override { return 3; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
private:
    std::vector<Index> nodes_;
};

// ═══ Face (2D) ═══
class Face : public Element {
public:
    int dimension() const override { return 2; }
    
    // 边界边 (按逆时针或右手定则)
    virtual std::vector<std::array<Index,2>> boundary_edges() const = 0;
};

// 3节点三角形
class Tri3 : public Face {
public:
    Tri3(Index n0, Index n1, Index n2) : nodes_{n0, n1, n2} {}
    
    ElementType type() const override { return ElementType::Tri3; }
    int num_nodes() const override { return 3; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 6节点三角形
class Tri6 : public Face {
public:
    Tri6(const std::array<Index,6>& ns) 
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Tri6; }
    int num_nodes() const override { return 6; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2, n3,n4,n5]
};

// 4节点四边形
class Quad4 : public Face {
public:
    Quad4(Index n0, Index n1, Index n2, Index n3) 
        : nodes_{n0, n1, n2, n3} {}
    
    ElementType type() const override { return ElementType::Quad4; }
    int num_nodes() const override { return 4; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[3]},
                {nodes_[3], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 8节点四边形
class Quad8 : public Face {
public:
    Quad8(const std::array<Index,8>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Quad8; }
    int num_nodes() const override { return 8; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::array<Index,2>> boundary_edges() const override {
        return {{nodes_[0], nodes_[1]},
                {nodes_[1], nodes_[2]},
                {nodes_[2], nodes_[3]},
                {nodes_[3], nodes_[0]}};
    }
    
private:
    std::vector<Index> nodes_;  // [n0,n1,n2,n3, n4,n5,n6,n7]
};

// ═══ Volume (3D) ═══
class Volume : public Element {
public:
    int dimension() const override { return 3; }
    
    // 边界面 (每个面由节点列表定义)
    virtual std::vector<std::vector<Index>> boundary_faces() const = 0;
};

// 4节点四面体
class Tet4 : public Volume {
public:
    Tet4(Index n0, Index n1, Index n2, Index n3)
        : nodes_{n0, n1, n2, n3} {}
    
    ElementType type() const override { return ElementType::Tet4; }
    int num_nodes() const override { return 4; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[2], nodes_[1]},
                {nodes_[0], nodes_[1], nodes_[3]},
                {nodes_[1], nodes_[2], nodes_[3]},
                {nodes_[2], nodes_[0], nodes_[3]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 10节点四面体
class Tet10 : public Volume {
public:
    Tet10(const std::array<Index,10>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Tet10; }
    int num_nodes() const override { return 10; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[2], nodes_[1]},
                {nodes_[0], nodes_[1], nodes_[3]},
                {nodes_[1], nodes_[2], nodes_[3]},
                {nodes_[2], nodes_[0], nodes_[3]}};
    }
    
private:
    std::vector<Index> nodes_;  // 4 corner + 6 edge mid
};

// 8节点六面体
class Brick8 : public Volume {
public:
    Brick8(const std::array<Index,8>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Brick8; }
    int num_nodes() const override { return 8; }
    int order() const override { return 1; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[3], nodes_[2], nodes_[1]},
                {nodes_[4], nodes_[5], nodes_[6], nodes_[7]},
                {nodes_[0], nodes_[1], nodes_[5], nodes_[4]},
                {nodes_[2], nodes_[3], nodes_[7], nodes_[6]},
                {nodes_[0], nodes_[4], nodes_[7], nodes_[3]},
                {nodes_[1], nodes_[2], nodes_[6], nodes_[5]}};
    }
    
private:
    std::vector<Index> nodes_;
};

// 20节点六面体
class Brick20 : public Volume {
public:
    Brick20(const std::array<Index,20>& ns)
        : nodes_(ns.begin(), ns.end()) {}
    
    ElementType type() const override { return ElementType::Brick20; }
    int num_nodes() const override { return 20; }
    int order() const override { return 2; }
    
    const std::vector<Index>& nodes() const override { return nodes_; }
    
    std::vector<std::vector<Index>> boundary_faces() const override {
        return {{nodes_[0], nodes_[3], nodes_[2], nodes_[1]},
                {nodes_[4], nodes_[5], nodes_[6], nodes_[7]},
                {nodes_[0], nodes_[1], nodes_[5], nodes_[4]},
                {nodes_[2], nodes_[3], nodes_[7], nodes_[6]},
                {nodes_[0], nodes_[4], nodes_[7], nodes_[3]},
                {nodes_[1], nodes_[2], nodes_[6], nodes_[5]}};
    }
    
private:
    std::vector<Index> nodes_;  // 8 corner + 12 edge mid
};

}  // namespace v2
}  // namespace fem
