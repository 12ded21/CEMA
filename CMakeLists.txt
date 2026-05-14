cmake_minimum_required(VERSION 3.12)  # 需支持C++17
project(CEMA LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CPLEX is an external dependency. Configure it with either:
#   cmake -S . -B build -DCPLEX_ROOT=/path/to/CPLEX2211
# or:
#   export CPLEX_ROOT=/path/to/CPLEX2211
set(CPLEX_ROOT "" CACHE PATH "Root directory of IBM ILOG CPLEX Optimization Studio")
if(NOT CPLEX_ROOT AND DEFINED ENV{CPLEX_ROOT})
    set(CPLEX_ROOT "$ENV{CPLEX_ROOT}" CACHE PATH "Root directory of IBM ILOG CPLEX Optimization Studio" FORCE)
endif()

if(NOT CPLEX_ROOT)
    message(FATAL_ERROR "CPLEX_ROOT is not set. Use -DCPLEX_ROOT=/path/to/CPLEX2211 or export CPLEX_ROOT=/path/to/CPLEX2211.")
endif()

set(CPLEX_INCLUDE_DIRS 
    ${CPLEX_ROOT}/cplex/include
    ${CPLEX_ROOT}/concert/include)
set(CPLEX_LIB_DIRS
    ${CPLEX_ROOT}/cplex/lib/x86-64_linux/static_pic
    ${CPLEX_ROOT}/concert/lib/x86-64_linux/static_pic)

if(NOT EXISTS "${CPLEX_ROOT}/cplex/include/ilcplex/ilocplex.h")
    message(FATAL_ERROR "Cannot find CPLEX headers under CPLEX_ROOT: ${CPLEX_ROOT}")
endif()

# 添加可执行文件
add_executable(main
    main.cpp
    ilp.cpp
    utils.cpp
    cliques.cpp
    get_mxclq.cpp)

# 包含头文件目录
target_include_directories(main PRIVATE ${CPLEX_INCLUDE_DIRS})

# 链接库目录
target_link_directories(main PRIVATE ${CPLEX_LIB_DIRS})

# 添加编译选项 (编译器插桩)
target_compile_options(main PRIVATE
    -O3
    -m64
    -fPIC
    -fno-strict-aliasing
    -fexceptions
    -w
    -DMOMC
    )  # 禁用警告（与-w等效）

# 添加预处理器定义
target_compile_definitions(main PRIVATE
    NDEBUG
    IL_STD)

# 链接库 (修复点：在这里加上链接器需要的 sanitizer 标志)
target_link_libraries(main PRIVATE
    concert
    ilocplex
    cplex
    pthread
    m
)
