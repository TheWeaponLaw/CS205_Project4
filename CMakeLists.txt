cmake_minimum_required(VERSION 3.16)

project(CS205_Project4)



aux_source_directory(./src DIR_SRCS)

include_directories(include)

add_executable(mat ${DIR_SRCS})

