# ── 示例程序 ──
add_executable(poisson_2d_v2 poisson_2d_v2.cpp)
target_link_libraries(poisson_2d_v2 PRIVATE fem)
target_include_directories(poisson_2d_v2 PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(heat_2d heat_2d.cpp)
target_link_libraries(heat_2d PRIVATE fem)
target_include_directories(heat_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_2d elasticity_2d.cpp)
target_link_libraries(elasticity_2d PRIVATE fem)
target_include_directories(elasticity_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_simple elasticity_simple.cpp)
target_link_libraries(elasticity_simple PRIVATE fem)
target_include_directories(elasticity_simple PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_debug elasticity_debug.cpp)
target_link_libraries(elasticity_debug PRIVATE fem)
target_include_directories(elasticity_debug PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_full_debug elasticity_full_debug.cpp)
target_link_libraries(elasticity_full_debug PRIVATE fem)
target_include_directories(elasticity_full_debug PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(cantilever_beam cantilever_beam.cpp)
target_link_libraries(cantilever_beam PRIVATE fem)
target_include_directories(cantilever_beam PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(thermal_stress_2d thermal_stress_2d.cpp)
target_link_libraries(thermal_stress_2d PRIVATE fem)
target_include_directories(thermal_stress_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)