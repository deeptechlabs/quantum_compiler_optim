cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

#Bring the headers, such as Student.h into the project
include_directories(include)
 
#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
 
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")

SET(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -D_REENTRANT -DVERSION=\"$(VERSION)\" -std=c++11")

add_executable(ibm_qx_mapping ${SOURCES})
