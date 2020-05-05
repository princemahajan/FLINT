#################################################################################
# Compiler Flags Module
# Author: Bharat Mahajan
#################################################################################

# Compiler flag check modules
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)

#[[#################################################################################

GetPlatformCompilerFlag

This function takes a list of compiler flags as an input, checks each of them against
the language compiler specified by the user and sets _Flag with the
correct flag that works, only if the current build configuration matches the provided
value. If the required flag is set to true, then it raises an 
error if the check fails otherwise it unsets the _Flag. 

Inputs:
Flag:           List of flags
IsRequired:     True if the flag must be set
LANG:           Language compiler to be checked against
Config:         The build Config in which the flag is used. Only 2 values
                are supported: Debug and Release.

Outputs:
_Flag:          The correct flag if found

One usage example:

set(COption /Zi -g)
GetPlatformCompilerFlag("${COption}" True Fortran)
target_compile_options(<target> PUBLIC "${_Flag}")

#################################################################################]]

function (GetPlatformCompilerFlag Flags IsRequired LANG)

# check each option one by one and return the first one that works
foreach(flag ${Flags})

    message(STATUS "GetPlatformCompilerFlag: checking option: " ${flag})

    # unset the variable in cache from the previous runs
    unset(FLAG_WORKS CACHE)

    # check the option against each language compiler
    if(LANG STREQUAL C)
        check_c_compiler_flag("${flag}" FLAG_WORKS)
    elseif(LANG STREQUAL CXX)
        check_cxx_compiler_flag("${flag}" FLAG_WORKS)
    elseif(LANG STREQUAL Fortran)
        check_fortran_compiler_flag("${flag}" FLAG_WORKS)
    else()
        message(FATAL_ERROR "GetPlatformCompilerFlag: Unknown language " ${LANG})
        break()
    endif()
    
    # if option works, set output and stop testing more flags
    if(FLAG_WORKS)
        set(_Flag ${flag} PARENT_SCOPE)
        set(FLAG_FOUND TRUE)        
        # unset the variables
        unset(FLAG_WORKS CACHE)
        break()
    else()
        set(FLAG_FOUND FALSE)
        message("GetPlatformCompilerFlag: Unknown flag: " ${flag})
        unset(_Flag PARENT_SCOPE)
    endif()
    
endforeach()

# if option is required and not found in list, raise an error
if(IsRequired AND (NOT FLAG_FOUND))
    message(FATAL_ERROR "GetPlatformCompilerFlag: Wrong required ${LANG} compiler flag " ${flag})
elseif((NOT IsRequired) AND (NOT FLAG_FOUND))
    # flag doesnt work but its ok!
    message(STATUS "GetPlatformCompilerFlag: Wrong ${LANG} compiler flag " ${flag})
endif()

endfunction()

