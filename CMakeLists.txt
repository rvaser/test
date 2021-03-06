cmake_minimum_required(VERSION 3.11)

project(test VERSION 0.0.4
             LANGUAGES CXX
             DESCRIPTION "Test.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

find_package(ram 2.0.0 QUIET)
if (NOT ram_FOUND)
  if (racon_install)
    message(FATAL_ERROR "Missing package ram!")
  endif ()

  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/lbcb-sci/ram
    GIT_TAG 3045360)

  FetchContent_GetProperties(ram)
  if (NOT ram_POPULATED)
    FetchContent_Populate(ram)
    add_subdirectory(
      ${ram_SOURCE_DIR}
      ${ram_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.13)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

add_executable(test
  src/main.cpp)

target_link_libraries(test
  bioparser::bioparser
  ram::ram)
