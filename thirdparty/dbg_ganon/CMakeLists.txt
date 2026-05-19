# =============================================================================
# ganon
# =============================================================================

# ！！！把 3.4 改成 >= 3.20，避免 “compatibility with CMake < 3.5 removed” 报错
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ganon VERSION 2.1.0 LANGUAGES CXX)

# ----------------------------------------------------------------------------- 
# build setup
# -----------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(FATAL_ERROR "Only GNU compilers are supported.")
endif()

option(VERBOSE_CONFIG "Verbose mode for quick build setup debugging" OFF)
option(CONDA "Flag for compilation in conda env." ON)         # ★ 在 conda 中构建，避免 -static -march=native
option(LONGREADS "Use uint32_t for very long reads" OFF)
option(INCLUDE_DIRS "Extra include dirs" "")

# ----------------------------------------------------------------------------- 
# build types
# -----------------------------------------------------------------------------

get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(isMultiConfig)
    if(NOT "Coverage" IN_LIST CMAKE_CONFIGURATION_TYPES)
        list(APPEND CMAKE_CONFIGURATION_TYPES Coverage)
    endif()
else()
    set(allowableBuildTypes Debug Release RelWithDebInfo MinSizeRel Coverage)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${allowableBuildTypes}")
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    elseif(NOT CMAKE_BUILD_TYPE IN_LIST allowableBuildTypes)
        message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}")
    endif()
endif()

# 覆盖 O3 等由外层注入的 flags 即可，不强制写死，避免和 conda 的 C(XX)FLAGS 冲突
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Coverage（保留即可）
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} -O0 --coverage")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -lgcov --coverage")

# warnings
add_compile_options(-Wall -Wextra -Wshadow -Wuninitialized -Wcast-align -Wunused
                    -Woverloaded-virtual -Wpedantic -Wnull-dereference -Wdouble-promotion
                    -Wformat=2 -Wstrict-aliasing -Wunused-variable -Wno-interference-size -Wno-null-dereference)
add_compile_options(-Wno-shadow -Wno-old-style-cast)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast)
endif()

# ★ 在 conda 环境中构建，不加 -static -march=native（原仓库里是 if(NOT CONDA) 才加）
if(LONGREADS)
    add_compile_options(-DLONGREADS)
endif()

# ----------------------------------------------------------------------------- 
# dependencies and 3rd party libraries
# -----------------------------------------------------------------------------

find_package(Threads REQUIRED)

# cxxopts（用 vendor）
add_library(cxxopts INTERFACE)
target_include_directories(cxxopts SYSTEM INTERFACE libs/cxxopts/include)

find_package(BZip2 REQUIRED)
find_package(ZLIB REQUIRED)

# ★ 关键：强制使用 vendor 的 seqan3（它自带 build_system）
#    不再依赖 conda 的 seqan3；通过在 CMAKE_PREFIX_PATH 里插入 vendor 的 build_system 来让 find_package 命中
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs/seqan3/build_system")
# 如果你的 GCC 比 seqan3 官方支持新（例如 13+），建议加此定义，抑制诊断
add_compile_definitions(SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC=1)
find_package(seqan3 3.3.0 REQUIRED)

# robin-hood-hashing（用 vendor）
add_library(robin-hood-hashing INTERFACE)
target_include_directories(robin-hood-hashing SYSTEM INTERFACE libs/robin-hood-hashing/src/include)

# ----------------------------------------------------------------------------- 
# verbose (可选)
# -----------------------------------------------------------------------------

if(VERBOSE_CONFIG)
    message(STATUS "  CMAKE_BUILD_TYPE        : ${CMAKE_BUILD_TYPE}")
    message(STATUS "  CMAKE_INSTALL_PREFIX    : ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "  CMAKE_PREFIX_PATH       : ${CMAKE_PREFIX_PATH}")
    message(STATUS "  SEQAN3_VERSION          : ${SEQAN3_VERSION}")
endif()

# ----------------------------------------------------------------------------- 
# targets
# -----------------------------------------------------------------------------

# ★ 完全禁用 tests：不要 enable_testing()/add_subdirectory(tests)
add_subdirectory(src)
