cmake_minimum_required(VERSION 3.16)
project(ggbp LANGUAGES CXX)

option(CXX "enable C++ compilation" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# define onde ficam os executáveis e bibliotecas
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})


option(BUILD_TESTS "Build unit tests" ON)

find_package(GUROBI REQUIRED)
find_package(Boost REQUIRED CONFIG)
find_package(OpenMP REQUIRED)

add_subdirectory(src)
add_subdirectory(app)

if(BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.11.0
  )
  FetchContent_MakeAvailable(catch2)

  add_subdirectory(tests)
endif()