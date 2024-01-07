# ____| |    _ _|  \  |__ __| 
# |     |      |    \ |   |   
# __|   |      |  |\  |   |   
# _|    _____|___|_| \_|  _|   
                                                        
# Author: Bharat Mahajan

cmake_minimum_required(VERSION 3.12)

###################
# FLINT Build Setup 
###################

project(FLINT VERSION 0.9 LANGUAGES Fortran)
enable_language(Fortran)

# set major and minor version
set(FLINT_MAJOR_VERSION 0)
set(FLINT_MINOR_VERSION 9)
set(FLINT_VERSION  ${FLINT_MAJOR_VERSION}.${FLINT_MINOR_VERSION})

# CMake module directory
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Prompt for user options

# Ask for the directory to install FLINT lib
set(FLINT_INSTALL_LIB_DIR lib/FLINT CACHE PATH "Installation directory for FLINT library")

# Ask for the directory to install TestFLINT executable
set(FLINT_INSTALL_BIN_DIR bin CACHE PATH "Installation directory for FLINT test executable")

# Make relative paths absolute
foreach(p LIB BIN)
  set(var FLINT_INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Lib dir should contain all the .mod files 
set(CMAKE_Fortran_MODULE_DIRECTORY ${FLINT_INSTALL_LIB_DIR})

# compile time constants (provide -D<VAR>=TRUE to cmake)
if (_REAL32)
  add_compile_definitions(_REAL32)
elseif (_REAL128)
  add_compile_definitions(_REAL128)
endif()

###############
# FLINT Targets
###############

set(FLINTBIN FLINT)

# Define the static library target
file(GLOB FLINT_SRC src/*.f90 src/*.F90)
add_library(${FLINTBIN} STATIC ${FLINT_SRC})
target_include_directories(${FLINTBIN} INTERFACE ${FLINT_INSTALL_LIB_DIR})

# The test executable for Work precision tests 
set(TestFLINT_WP_SRC tests/test_WP.f90)
set(TestFLINT_WPBIN TestFLINT_WP)
add_executable(${TestFLINT_WPBIN} ${TestFLINT_WP_SRC})

# link
target_link_libraries(${TestFLINT_WPBIN} ${FLINTBIN})

# The test executable for CR3BP test
set(TestFLINT_SRC tests/test_CR3BP.f90)
set(TestFLINTBIN TestFLINT_CR3BP)
add_executable(${TestFLINTBIN} ${TestFLINT_SRC})

# link
target_link_libraries(${TestFLINTBIN} ${FLINTBIN})

# The test executable for Lorenz diffeq
set(TestFLINT_Lorenz_SRC tests/test_Lorenz.f90)
set(TestFLINT_LorenzBIN TestFLINT_Lorenz)
add_executable(${TestFLINT_LorenzBIN} ${TestFLINT_Lorenz_SRC})

# link
target_link_libraries(${TestFLINT_LorenzBIN} ${FLINTBIN})

# The test executable for TBP test
set(TestFLINT_SRC tests/test_TBP.f90)
set(TestFLINTBIN TestFLINT_TBP)
add_executable(${TestFLINTBIN} ${TestFLINT_SRC})

# link
target_link_libraries(${TestFLINTBIN} ${FLINTBIN})

###########################
# Compile Options for FLINT
###########################

include(CompilerOptions)

# Thee options are for gfortran
set(COption "-std=f2008ts")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
add_compile_options("${_Flag}")
set(COption "-ffree-line-length-0")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
add_compile_options("${_Flag}")

# Add ifort specific options

# standard
set(COption "/stand:f08" "-stand f08")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PUBLIC "${_Flag}")

# Data alignment always
set(COption "/align:all" "-align all")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PUBLIC "${_Flag}")

# Traceback information
set(COption "/traceback" "-traceback")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE ${_Flag})

# Generate auto-dispatch code paths for better optimization if available 
set(COption "/QxHost" "-xHost")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE ${_Flag})

# # optimization suggestions from compiler (doesn't create excecutables)
# set(COption "/Qguide" "-guide")
# GetPlatformCompilerFlag("${COption}" FALSE Fortran)
# target_compile_options(${FLINTBIN} PRIVATE ${_Flag}>)

# enable changing of the memory lauyut
set(COption "/Qpad" "-pad")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE ${_Flag})

set(COption "/warn:all" "-warn all")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE ${_Flag})

#######################
# Debug Compile Options
#######################

# No optimization
set(COption "/Od" "-O0")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Debug>:${_Flag}>)

# Generate full debug information, passed to the linker
set(COption "/debug:full" "-debug full")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Debug>:${_Flag}>)

# Generate full debug information in object or PDB file, passed to compiler
set(COption "/Zi" "-g2")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Debug>:${_Flag}>)

# Debug information for Parameters
set(COption "/debug-parameters" "-debug-parameters")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Debug>:${_Flag}>)

# Array bounds checking
set(COption "/check:all" "-check all" "-fcheck=all")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Debug>:${_Flag}>)


#########################
# Release Compile Options
#########################

# Aggressive optimization for maximum speed: includes bunch of options
# Make sure we also set the appropriat Intel instruction set to use after this
set(COption "/fast" "-fast")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Release>:${_Flag}>)

# Most optimization enabled
set(COption "/O3" "-O3")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Release>:${_Flag}>)

# Enable auto-parallelizer
set(COption "/Qparallel" "-parallel")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Release>:${_Flag}>)

# Enable multi-file optimization
set(COption "/Qipo" "-ipo")
GetPlatformCompilerFlag("${COption}" FALSE Fortran)
target_compile_options(${FLINTBIN} PRIVATE $<$<CONFIG:Release>:${_Flag}>)


############################
# Install and Export Targets
############################


install(TARGETS ${FLINTBIN}
  # IMPORTANT: Add the FLINT library to the "export-set"
  EXPORT FLINTTargets
  LIBRARY DESTINATION "${FLINT_INSTALL_LIB_DIR}" COMPONENT lib
  ARCHIVE DESTINATION "${FLINT_INSTALL_LIB_DIR}" COMPONENT lib
  COMPONENT dev)

install(TARGETS ${TestFLINTBIN}
  # IMPORTANT: Add the TestFLINT executable to the "export-set"
  EXPORT FLINTTargets
  RUNTIME DESTINATION "${FLINT_INSTALL_BIN_DIR}" COMPONENT bin)

# Build-tree export set
export(EXPORT FLINTTargets 
    NAMESPACE FLINT::
    FILE "${FLINT_INSTALL_LIB_DIR}/cmake/FLINTTargets.cmake")

##########################
# Config and Version Files
##########################

include(CMakePackageConfigHelpers)

# version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/FLINTConfigVersion.cmake
    VERSION ${FLINT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FLINTConfig.cmake.in 
        ${CMAKE_CURRENT_BINARY_DIR}/FLINTConfig.cmake
        INSTALL_DESTINATION cmake
        INSTALL_PREFIX ${FLINT_INSTALL_LIB_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FLINTConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FLINTConfigVersion.cmake
        DESTINATION ${FLINT_INSTALL_LIB_DIR}/cmake)