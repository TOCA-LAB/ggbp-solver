set(GGBP_SOURCES
    graph.cpp
    ggbp_utils.cpp
    greedy.cpp
    grasp.cpp
    cover.cpp
    naive.cpp
    prym.cpp
    pap.cpp
)

add_library(ggbp STATIC ${GGBP_SOURCES}
    brkga.cpp)
target_include_directories(ggbp
    PUBLIC ${PROJECT_SOURCE_DIR}/include ${GUROBI_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
    PRIVATE ${PROJECT_SOURCE_DIR}/src
)
add_library(gurobi STATIC IMPORTED)
if(MSVC)
    target_compile_options(ggbp PRIVATE -W1)
else()
    target_compile_options(ggbp PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(ggbp
    optimized ${GUROBI_CXX_LIBRARY}
    debug ${GUROBI_CXX_DEBUG_LIBRARY}
)

target_link_libraries(ggbp ${GUROBI_LIBRARY})