cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project("tmath")
set(CMAKE_CXX_STANDARD 20)

add_compile_options(
	-Wall
	-Wextra
	-Wno-missing-braces
	-O3
)

add_executable(testexc
test.cpp
)
