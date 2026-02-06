# ── 示例程序 ──
add_executable(poisson_2d poisson_2d.cpp)
target_link_libraries(poisson_2d PRIVATE fem)
target_include_directories(poisson_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(heat_conduction_2d heat_conduction_2d.cpp)
target_link_libraries(heat_conduction_2d PRIVATE fem)
target_include_directories(heat_conduction_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(elasticity_2d elasticity_2d.cpp)
target_link_libraries(elasticity_2d PRIVATE fem)
target_include_directories(elasticity_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(thermal_stress_2d thermal_stress_2d.cpp)
target_link_libraries(thermal_stress_2d PRIVATE fem)
target_include_directories(thermal_stress_2d PRIVATE ${CMAKE_SOURCE_DIR}/src)
