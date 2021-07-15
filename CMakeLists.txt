cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(fem VERSION 1.0)


enable_testing()

include_directories(${PROJECT_SOURCE_DIR} args-parser)
add_executable(fem main.cpp interpol.cpp utils.cpp builders.cpp)

target_link_libraries(fem gtest)
target_link_libraries(fem pthread)
target_link_libraries(fem lapacke)
target_link_libraries(fem lapack)
target_link_libraries(fem cblas)
target_link_libraries(fem blas)


include(GoogleTest)
gtest_discover_tests(fem)