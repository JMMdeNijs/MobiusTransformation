cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project("MobiusTransformation")

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -pthread")
find_package (Threads)

find_package(X11 REQUIRED)
message(STATUS "X11_FOUND = ${X11_FOUND}")
message(STATUS "X11_INCLUDE_DIR = ${X11_INCLUDE_DIR}")
message(STATUS "X11_LIBRARIES = ${X11_LIBRARIES}")

include_directories(	
			~/Libraries/CImg
			~/Libraries/Eigen32
			${X11_INCLUDE_DIR}	)

add_executable(MobiusTransformation 
			mobius_Transformation_v1.0.cpp
			mobiuslib_v1.0.cpp )

target_link_libraries (	MobiusTransformation
			${CMAKE_THREAD_LIBS_INIT}
			${X11_LIBRARIES})
