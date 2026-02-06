# ── GoogleTest ──
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

# ── 测试目标 (按模块逐步添加) ──
add_executable(fem_tests
    test_core.cpp
    test_mesh.cpp
    test_element.cpp
    test_assembly.cpp
    test_solver.cpp
    test_io.cpp
    test_physics.cpp
)

target_link_libraries(fem_tests PRIVATE fem gtest_main)
target_include_directories(fem_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME fem_tests COMMAND fem_tests)
