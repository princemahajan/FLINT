# ____| |    _ _|  \  |__ __| 
# |     |      |    \ |   |   
# __|   |      |  |\  |   |   
# _|    _____|___|_| \_|  _|   
                                                        

# Version: 0.9
# Author: Bharat Mahajan

cmake_minimum_required(VERSION 3.12)

# Project
project(FLINT VERSION 0.9 LANGUAGES Fortran)
enable_language(Fortran)

# Prompt for user options

# Extract Fortran Compiler name
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# TestFLINT executable
set(TestFLINTBIN TestFLINT)

# TestFLINT static library
set(FLINTBIN flint)

# Define the static library target
file(GLOB FLINT_SRC src/*.f90)
add_library(${FLINTBIN} STATIC ${FLINT_SRC})

# Define the executable
file(GLOB TestFLINT_SRC src/tests/*.f90)
add_executable(${TestFLINTBIN} ${TestFLINT_SRC})

# add gfortran specific options
if (${Fortran_COMPILER_NAME} STREQUAL "gfortran.exe" OR ${Fortran_COMPILER_NAME} STREQUAL "gfortran")
    target_compile_options(${FLINTBIN} PRIVATE "-std=f2008ts")
    target_compile_options(${FLINTBIN} PRIVATE "-ffree-line-length-0")
    target_compile_options(${TestFLINTBIN} PRIVATE "-std=f2008ts")
    target_compile_options(${TestFLINTBIN} PRIVATE "-ffree-line-length-0")    
endif()

# link flint
target_link_libraries(${TestFLINTBIN} ${FLINTBIN})


