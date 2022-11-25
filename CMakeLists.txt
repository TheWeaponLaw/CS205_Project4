cmake_minimum_required(VERSION 3.16)

project(CS205_Project4)

aux_source_directory(./src DIR_SRCS)

include_directories(include)# /usr/include/)

#link_libraries(-lopenblas -lpthread -lgfortran)

#target_link_libraries(mat /usr/local/OpenBLAS/lib/libopenblas.a)
#target_link_libraries(mat pthread)
#target_link_libraries(mat cblas)

add_executable(mat ${DIR_SRCS})

target_link_libraries(mat /usr/local/OpenBLAS/lib/libopenblas.a)
target_link_libraries(mat pthread)

add_definitions(-mavx2)

add_definitions(-fopenmp -lgomp)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
   # target_link_libraries(mat PUBLIC OpenMP::OpenMP_CXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#cmake . -DCMAKE_BUILD_TYPE=Release -Bbuild
#-I /usr/include/ -L/usr/lib -lopenblas -lpthread -lgfortran

