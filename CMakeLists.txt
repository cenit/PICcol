cmake_minimum_required (VERSION 2.8)
project (PICcol)
set(CMAKE_CXX_STANDARD 11)

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

include_directories(src/jsoncons/include)

add_executable(PICcol
  src/datatypes.cpp
  src/datatypes.h
  src/evolution.cpp
  src/evolution.h
  src/filler.cpp
  src/filler.h
  src/main.cpp
)

install(TARGETS PICcol DESTINATION bin)
