if(NOT dune-integrate_FOUND)
# Whether this module is installed or not
set(dune-integrate_INSTALLED ON)

# Settings specific to the module

# Package initialization

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was dune-integrate-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#report other information
set_and_check(dune-integrate_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-integrate_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set(dune-integrate_CXX_FLAGS "-std=c++17 ")
set(dune-integrate_CXX_FLAGS_DEBUG "-g")
set(dune-integrate_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-integrate_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-integrate_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-integrate_DEPENDS "dune-common;dune-geometry;dune-grid;dune-istl;dune-localfunctions;dune-functions")
set(dune-integrate_SUGGESTS "")
set(dune-integrate_MODULE_PATH "${PACKAGE_PREFIX_DIR}/share/dune/cmake/modules")
set(dune-integrate_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-integrate_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-integrate-targets.cmake")
endif()
endif()
