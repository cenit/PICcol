cmake_minimum_required (VERSION 2.8)
project (PICcol)
set(CMAKE_CXX_STANDARD 11)

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

include_directories(${CMAKE_SOURCE_DIR}/src/jsoncons/include)

add_executable(PICcol
  ${CMAKE_SOURCE_DIR}/src/datatypes.cpp
  ${CMAKE_SOURCE_DIR}/src/datatypes.h
  ${CMAKE_SOURCE_DIR}/src/evolution.cpp
  ${CMAKE_SOURCE_DIR}/src/evolution.h
  ${CMAKE_SOURCE_DIR}/src/filler.cpp
  ${CMAKE_SOURCE_DIR}/src/filler.h
  ${CMAKE_SOURCE_DIR}/src/main.cpp
)
