# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/munna/.local/lib/python3.8/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/munna/.local/lib/python3.8/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/munna/dune/dune-integrate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/munna/dune/dune-integrate/build-cmake

# Utility rule file for OUTPUT.

# Include the progress variables for this target.
include CMakeFiles/OUTPUT.dir/progress.make

CMakeFiles/OUTPUT:
	config_collected.h.cmake
	dune_regenerate_config_cmake ( )

OUTPUT: CMakeFiles/OUTPUT
OUTPUT: CMakeFiles/OUTPUT.dir/build.make

.PHONY : OUTPUT

# Rule to build all files generated by this target.
CMakeFiles/OUTPUT.dir/build: OUTPUT

.PHONY : CMakeFiles/OUTPUT.dir/build

CMakeFiles/OUTPUT.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/OUTPUT.dir/cmake_clean.cmake
.PHONY : CMakeFiles/OUTPUT.dir/clean

CMakeFiles/OUTPUT.dir/depend:
	cd /home/munna/dune/dune-integrate/build-cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/munna/dune/dune-integrate /home/munna/dune/dune-integrate /home/munna/dune/dune-integrate/build-cmake /home/munna/dune/dune-integrate/build-cmake /home/munna/dune/dune-integrate/build-cmake/CMakeFiles/OUTPUT.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/OUTPUT.dir/depend

