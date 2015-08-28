cmake_minimum_required(VERSION 3.2)
project(TheiaSolver)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib )

include_directories(
${PROJECT_SOURCE_DIR}/include


)

set(SOURCE_FILES )

# Add sources
set(THEIA_SRC
  math/closed_form_polynomial_solver.cc
  math/find_polynomial_roots_companion_matrix.cc
  math/find_polynomial_roots_jenkins_traub.cc
  math/matrix/dominant_eigensolver.cc
  math/polynomial.cc
  math/probability/sequential_probability_ratio.cc

  util/filesystem.cc
  util/random.cc
  util/stringprintf.cc
  util/threadpool.cc
  util/timer.cc
)

add_library( TheiaSolver ${THEIA_SRC} )
