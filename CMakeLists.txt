#!If your scip is in Release mode, set Release mode: cmake -DCMAKE_BUILD_TYPE=Release
# Otherwise your scip is in Debug mode, set Release mode: cmake -DCMAKE_BUILD_TYPE=Debug
cmake_minimum_required(VERSION 3.0)

# set build type
set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


project(wUMCFC)

# find SCIP libraries and include directories
find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})



add_executable(wUMCFC
src_bp/cppmain.cpp
src_bp/Graph.cpp
src_bp/ProblemGraph.cpp 
src_bp/PricingGraph.cpp 
src_bp/ProbDataNTRT.cpp 
src_bp/ReaderNTRT.cpp
src_bp/table_xyz.cpp
src_bp/ConshdlrNTRT.cpp
src_bp/BranchNTRT.cpp 
src_bp/PricerNTRT.cpp
src_bp/check.cpp
    )


# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()


target_link_libraries(wUMCFC -lm -lpthread -ldl ${SCIP_LIBRARIES} ${LIBM})
