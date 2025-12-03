cmake_minimum_required(VERSION 3.12...3.31)
project(MyTeamCGALProject)

# Point to CGAL submodule
set(CGAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cgal)

find_package(CGAL REQUIRED)

# Add your custom examples
add_subdirectory(examples)
