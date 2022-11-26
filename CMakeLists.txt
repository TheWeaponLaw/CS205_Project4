cmake_minimum_required(VERSION 3.16)

set(CMAKE_C_STANDARD 11)
project(CS205_Project4)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
aux_source_directory(./src DIR_SRCS)

include_directories(include)

add_executable(mat ${DIR_SRCS})

target_link_libraries(mat /usr/local/OpenBLAS/lib/libopenblas.a)
target_link_libraries(mat pthread)                        

#add_compile_options(-O3)

add_definitions(-mavx2)

add_definitions(-fopenmp -lgomp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
    target_link_libraries(mat OpenMP::OpenMP_CXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#cmake . -DCMAKE_BUILD_TYPE=Release -Bbuild

