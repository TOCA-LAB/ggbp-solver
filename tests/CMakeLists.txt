enable_testing()

include(Catch)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_ggbp
  test_graph.cpp
  test_combinations.cpp
  test_greedy.cpp
  test_regression_graph_io.cpp
  test_regression_graph_cache.cpp
  test_regression_graph_contracts.cpp
  test_regression_graph_isolated.cpp
  test_regression_graph_simulation.cpp
  test_regression_heuristics_core.cpp
  test_regression_heuristics_isolated.cpp
  test_regression_heuristics_build.cpp
)

target_include_directories(test_ggbp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_ggbp PRIVATE ggbp Catch2::Catch2WithMain)

catch_discover_tests(test_ggbp
  TEST_PREFIX "core."
  PROPERTIES
    LABELS "regression-core"
    TIMEOUT 10
)

add_executable(test_ggbp_gurobi
  test_regression_models.cpp
)

target_include_directories(test_ggbp_gurobi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_ggbp_gurobi PRIVATE ggbp Catch2::Catch2WithMain)

catch_discover_tests(test_ggbp_gurobi
  TEST_PREFIX "gurobi."
  PROPERTIES
    LABELS "regression-gurobi"
    RESOURCE_LOCK gurobi
    TIMEOUT 60
)

# Findings whose safe regression is a build/source policy rather than an
# in-process behavior. Register one CTest per finding so failures remain
# independently actionable. These checks are red on the audited revision.
set(GGBP_SOURCE_REGRESSION_FINDINGS
  F01 F05 F06 F09 F10 F12 F14 F16 F18 F20 F21 F22 F24 F25 F26 F27
)

foreach(finding IN LISTS GGBP_SOURCE_REGRESSION_FINDINGS)
  add_test(
    NAME "source.${finding}"
    COMMAND
      ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/regression_source_checks.py
      ${finding}
      --source-root ${PROJECT_SOURCE_DIR}
  )
  set_tests_properties("source.${finding}" PROPERTIES
    LABELS "regression;source;${finding}"
    TIMEOUT 10
  )
endforeach()
