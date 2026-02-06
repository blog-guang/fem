# ── 示例程序 ──
add_executable(poisson_2d_demo poisson_2d_demo.cpp)
target_link_libraries(poisson_2d_demo PRIVATE fem)
target_include_directories(poisson_2d_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(poisson_2d_v2 poisson_2d_v2.cpp)
target_link_libraries(poisson_2d_v2 PRIVATE fem)
target_include_directories(poisson_2d_v2 PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(heat_2d heat_2d.cpp)
target_link_libraries(heat_2d PRIVATE fem)
target_include_directories(heat_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_2d elasticity_2d.cpp)
target_link_libraries(elasticity_2d PRIVATE fem)
target_include_directories(elasticity_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)