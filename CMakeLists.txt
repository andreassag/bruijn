cmake_minimum_required(VERSION 3.20)
project(bruijn LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_executable(bruijn_assembler
    src/main.cpp
    src/read_buffer.cpp
)

target_include_directories(bruijn_assembler PRIVATE include)

add_executable(bruijn_tests
    tests/test_graph.cpp
    src/read_buffer.cpp
)

target_include_directories(bruijn_tests PRIVATE include)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bruijn_assembler PRIVATE
        -O3
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
    target_compile_options(bruijn_tests PRIVATE
        -O3
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
endif()

add_test(NAME bruijn_tests COMMAND bruijn_tests)
