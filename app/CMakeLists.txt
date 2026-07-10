link_libraries(nlohmann_json::nlohmann_json argparse)

add_executable(cover cover.cpp)
target_link_libraries(cover PRIVATE ggbp)

add_executable(naive naive.cpp)
target_link_libraries(naive PRIVATE ggbp)

add_executable(prym prym.cpp)
target_link_libraries(prym PRIVATE ggbp)

add_executable(pap pap.cpp)
target_link_libraries(pap PRIVATE ggbp)

add_executable(brkga brkga.cpp )
target_link_libraries(brkga PRIVATE ggbp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(brkga PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(grasp grasp.cpp)
target_link_libraries(grasp PRIVATE ggbp)

add_executable(unburned_neighbors unburned_neighbors.cpp)
target_link_libraries(unburned_neighbors PRIVATE ggbp)

add_executable(greedy greedy.cpp)
target_link_libraries(greedy PRIVATE ggbp)

add_executable(dominating dominating.cpp)
target_link_libraries(dominating PRIVATE ggbp)

add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PRIVATE ggbp)
