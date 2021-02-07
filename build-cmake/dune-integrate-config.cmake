if(NOT dune-integrate_FOUND)
# Whether this module is installed or not
set(dune-integrate_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/munna/dune/dune-integrate)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-integrate_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-integrate_INCLUDE_DIRS "/home/munna/dune/dune-integrate")
set(dune-integrate_CXX_FLAGS "-std=c++17 ")
set(dune-integrate_CXX_FLAGS_DEBUG "-g")
set(dune-integrate_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-integrate_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-integrate_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-integrate_DEPENDS "dune-common;dune-geometry;dune-grid;dune-istl;dune-localfunctions;dune-functions")
set(dune-integrate_SUGGESTS "")
set(dune-integrate_MODULE_PATH "/home/munna/dune/dune-integrate/cmake/modules")
set(dune-integrate_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-integrate_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-integrate-targets.cmake")
endif()
endif()
