# ── 示例程序 ──
add_executable(poisson_2d_demo poisson_2d_demo.cpp)

target_link_libraries(poisson_2d_demo PRIVATE fem)
target_include_directories(poisson_2d_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)

# 如果需要更多示例，可以在这里添加
# add_executable(example2 example2.cpp)
# target_link_libraries(example2 PRIVATE fem)