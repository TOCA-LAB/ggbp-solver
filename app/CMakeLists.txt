add_executable(cover cover.cpp)
target_link_libraries(cover PRIVATE ggbp)

add_executable(naive naive.cpp)
target_link_libraries(naive PRIVATE ggbp)

add_executable(greedy greedy.cpp )
target_link_libraries(greedy PRIVATE ggbp)

add_executable(brkga brkga.cpp )
target_link_libraries(brkga PRIVATE ggbp)

add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PRIVATE ggbp)
