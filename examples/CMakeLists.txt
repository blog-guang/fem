# ── 示例程序 ──
add_executable(poisson_2d_demo poisson_2d_demo.cpp)
target_link_libraries(poisson_2d_demo PRIVATE fem)
target_include_directories(poisson_2d_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(poisson_2d_v2 poisson_2d_v2.cpp)
target_link_libraries(poisson_2d_v2 PRIVATE fem)
target_include_directories(poisson_2d_v2 PRIVATE ${CMAKE_SOURCE_DIR}/src)