cmake_minimum_required(VERSION 3.15)
project(MetaMaxwell)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fconstexpr-steps=1000000000 ${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "-fconstexpr-ops-limit=1000000000 ${CMAKE_CXX_FLAGS}")

add_executable(MetaMaxwell main.cpp)