# Install script for directory: /home/munna/dune/dune-integrate

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/munna/dune")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  set(CMAKE_MODULE_PATH /home/munna/dune/dune-vtk/cmake/modules;/home/munna/dune/dune-functions/cmake/modules;/home/munna/dune/dune-localfunctions/cmake/modules;/home/munna/dune/dune-grid/cmake/modules;/home/munna/dune/dune-istl/cmake/modules;/home/munna/dune/dune-typetree/cmake/modules;/home/munna/dune/dune-uggrid/cmake/modules;/home/munna/dune/dune-geometry/cmake/modules;/home/munna/dune/dune-common/cmake/modules;/home/munna/dune/dune-integrate/cmake/modules;/home/munna/dune/dune-common/cmake/modules/FindPkgConfig;/home/munna/dune/dune-common/cmake/modules/FindPython3)
              set(DUNE_PYTHON_WHEELHOUSE /home/munna/dune/share/dune/wheelhouse)
              include(DuneExecuteProcess)
              dune_execute_process(COMMAND "/home/munna/.local/lib/python3.8/site-packages/cmake/data/bin/cmake" --build . --target install_python --config $<CONFIG>)
              
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/dunecontrol/dune-integrate" TYPE FILE FILES "/home/munna/dune/dune-integrate/dune.module")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/dune-integrate" TYPE FILE FILES
    "/home/munna/dune/dune-integrate/build-cmake/cmake/pkg/dune-integrate-config.cmake"
    "/home/munna/dune/dune-integrate/build-cmake/dune-integrate-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/dune-integrate" TYPE FILE FILES "/home/munna/dune/dune-integrate/config.h.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/munna/dune/dune-integrate/build-cmake/dune-integrate.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/munna/dune/dune-integrate/build-cmake/src/cmake_install.cmake")
  include("/home/munna/dune/dune-integrate/build-cmake/dune/cmake_install.cmake")
  include("/home/munna/dune/dune-integrate/build-cmake/doc/cmake_install.cmake")
  include("/home/munna/dune/dune-integrate/build-cmake/cmake/modules/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/munna/dune/dune-integrate/build-cmake/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
