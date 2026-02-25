add_executable(cover cover.cpp)
target_link_libraries(cover PRIVATE ggbp)

add_executable(greedy greedy.cpp )
target_link_libraries(greedy PRIVATE ggbp)

add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PRIVATE ggbp)
