# Material module

# 源文件
set(MATERIAL_SOURCES
    state_variables.cpp
    material.cpp
)

# 头文件（用于IDE组织）
set(MATERIAL_HEADERS
    state_variables.h
    material.h
)

# 添加到父目标（fem库）
target_sources(fem PRIVATE ${MATERIAL_SOURCES})

# 安装头文件
install(FILES ${MATERIAL_HEADERS}
        DESTINATION include/fem/material)
