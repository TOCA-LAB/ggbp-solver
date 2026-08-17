set(GGBP_SOURCES
    graph.cpp
    ggbp_utils.cpp
    burning_simulation.cpp
    gurobi_status.cpp
    cover.cpp
    naive.cpp
    prym.cpp
    pap.cpp
    grasp.cpp
    brkga.cpp
    spontaneous_neighbors.cpp
    dominating.cpp
    unburned_neighbors.cpp
)

add_library(ggbp STATIC ${GGBP_SOURCES})

if(NOT ENABLE_CONTRACTS)
    target_compile_definitions(ggbp PUBLIC NO_CONTRACTS=1)
endif()

target_include_directories(ggbp
    PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS}
    PRIVATE ${PROJECT_SOURCE_DIR}/src
)

if(MSVC)
    target_compile_options(ggbp PRIVATE -W1)
else()
    target_compile_options(ggbp PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(ggbp
    PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
    PUBLIC
        OpenMP::OpenMP_CXX
        gurobi
)
