cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-g")
project(NTRT)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(ntrt_bb
   src_bb/cppmain.cpp
   src_bb/Graph.cpp
   src_bb/ProblemGraph.cpp 
   src_bb/ProbDataNTRT.cpp 
   src_bb/ReaderNTRT.cpp
   src_bb/BranchNTRT.cpp 
   src_bb/check.cpp
    )

target_link_libraries(ntrt_bb ${SCIP_LIBRARIES} ${LIBM})
