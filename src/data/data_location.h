/**
 * data_location.h - 数据存储位置枚举
 * 
 * 定义场数据在网格中的存储位置
 */

#pragma once

#include <string>

namespace fem {
namespace data {

/**
 * 数据存储位置
 */
enum class DataLocation {
    Node,        ///< 节点（用于连续场，如位移、温度）
    Element,     ///< 单元中心（用于单元平均值，如应力、应变）
    Face,        ///< 面（用于通量）
    Edge,        ///< 边（用于 2D 问题）
    GaussPoint,  ///< 高斯积分点（用于积分点数据）
    Unknown      ///< 未知位置
};

/**
 * 将 DataLocation 转换为字符串
 */
inline std::string to_string(DataLocation loc) {
    switch (loc) {
        case DataLocation::Node:       return "Node";
        case DataLocation::Element:    return "Element";
        case DataLocation::Face:       return "Face";
        case DataLocation::Edge:       return "Edge";
        case DataLocation::GaussPoint: return "GaussPoint";
        default:                       return "Unknown";
    }
}

} // namespace data
} // namespace fem
