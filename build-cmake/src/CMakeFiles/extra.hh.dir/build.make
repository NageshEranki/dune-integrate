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

# Include any dependencies generated for this target.
include src/CMakeFiles/extra.hh.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/extra.hh.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/extra.hh.dir/flags.make

src/CMakeFiles/extra.hh.dir/oned.cc.o: src/CMakeFiles/extra.hh.dir/flags.make
src/CMakeFiles/extra.hh.dir/oned.cc.o: ../src/oned.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munna/dune/dune-integrate/build-cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/extra.hh.dir/oned.cc.o"
	cd /home/munna/dune/dune-integrate/build-cmake/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/extra.hh.dir/oned.cc.o -c /home/munna/dune/dune-integrate/src/oned.cc

src/CMakeFiles/extra.hh.dir/oned.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/extra.hh.dir/oned.cc.i"
	cd /home/munna/dune/dune-integrate/build-cmake/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munna/dune/dune-integrate/src/oned.cc > CMakeFiles/extra.hh.dir/oned.cc.i

src/CMakeFiles/extra.hh.dir/oned.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/extra.hh.dir/oned.cc.s"
	cd /home/munna/dune/dune-integrate/build-cmake/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munna/dune/dune-integrate/src/oned.cc -o CMakeFiles/extra.hh.dir/oned.cc.s

# Object files for target extra.hh
extra_hh_OBJECTS = \
"CMakeFiles/extra.hh.dir/oned.cc.o"

# External object files for target extra.hh
extra_hh_EXTERNAL_OBJECTS =

src/libextra.hh.a: src/CMakeFiles/extra.hh.dir/oned.cc.o
src/libextra.hh.a: src/CMakeFiles/extra.hh.dir/build.make
src/libextra.hh.a: src/CMakeFiles/extra.hh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/munna/dune/dune-integrate/build-cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libextra.hh.a"
	cd /home/munna/dune/dune-integrate/build-cmake/src && $(CMAKE_COMMAND) -P CMakeFiles/extra.hh.dir/cmake_clean_target.cmake
	cd /home/munna/dune/dune-integrate/build-cmake/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/extra.hh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/extra.hh.dir/build: src/libextra.hh.a

.PHONY : src/CMakeFiles/extra.hh.dir/build

src/CMakeFiles/extra.hh.dir/clean:
	cd /home/munna/dune/dune-integrate/build-cmake/src && $(CMAKE_COMMAND) -P CMakeFiles/extra.hh.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/extra.hh.dir/clean

src/CMakeFiles/extra.hh.dir/depend:
	cd /home/munna/dune/dune-integrate/build-cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/munna/dune/dune-integrate /home/munna/dune/dune-integrate/src /home/munna/dune/dune-integrate/build-cmake /home/munna/dune/dune-integrate/build-cmake/src /home/munna/dune/dune-integrate/build-cmake/src/CMakeFiles/extra.hh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/extra.hh.dir/depend

