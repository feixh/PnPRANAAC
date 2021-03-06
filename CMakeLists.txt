cmake_minimum_required(VERSION 2.8)
project(TheiaSolver)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
set(CMAKE_MODULE_PATH $CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMakeModules )
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib )

find_package(Eigen3 REQUIRED)

find_package(Glog REQUIRED)

include_directories(
${PROJECT_SOURCE_DIR}/include
${EIGEN3_INCLUDE_DIR}
${GLOG_INCLUDE_DIRS}
)

# Add sources
set(RANSAC_SRC
  src/math/closed_form_polynomial_solver.cc
  src/math/find_polynomial_roots_companion_matrix.cc
  src/math/find_polynomial_roots_jenkins_traub.cc
  src/math/matrix/dominant_eigensolver.cc
  src/math/polynomial.cc
  src/math/probability/sequential_probability_ratio.cc
  src/util/random.cc
#   src/util/stringprintf.cc
#   src/util/threadpool.cc
  src/util/timer.cc

  src/pnpsolvers/P3P_Kneip.cpp
)

add_library( RANSAC STATIC ${RANSAC_SRC} )

link_libraries(
RANSAC
${GLOG_LIBRARIES}
)

add_executable( ransac_test test/ransac_test.cpp)

add_executable( p3p_test test/p3p_test.cpp)
