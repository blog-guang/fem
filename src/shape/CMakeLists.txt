# Shape function module

# 源文件
set(SHAPE_SOURCES
    shape_function.cpp
    shape_tri3.cpp
    shape_quad4.cpp
    shape_tet4.cpp
    shape_brick8.cpp
    shape_function_factory.cpp
)

# 头文件（用于IDE组织）
set(SHAPE_HEADERS
    shape_function.h
    shape_tri3.h
    shape_quad4.h
    shape_tet4.h
    shape_brick8.h
    shape_function_factory.h
)

# 添加到父目标（fem库）
target_sources(fem PRIVATE ${SHAPE_SOURCES})

# 安装头文件
install(FILES ${SHAPE_HEADERS}
        DESTINATION include/fem/shape)
