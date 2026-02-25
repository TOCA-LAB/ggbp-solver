enable_testing()

add_executable(test_ggbp
  test_graph.cpp
  test_combinations.cpp
  test_greedy.cpp
)

target_link_libraries(test_ggbp
  PRIVATE
    ggbp
    Catch2::Catch2WithMain
)

add_test(NAME graph_tests COMMAND test_ggbp)



