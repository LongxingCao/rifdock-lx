# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/longxing/Rifdock/rifdock

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/longxing/Rifdock/rifdock/build

# Include any dependencies generated for this target.
include external/gmock/gtest/CMakeFiles/gtest.dir/depend.make

# Include the progress variables for this target.
include external/gmock/gtest/CMakeFiles/gtest.dir/progress.make

# Include the compile flags for this target's objects.
include external/gmock/gtest/CMakeFiles/gtest.dir/flags.make

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o: external/gmock/gtest/CMakeFiles/gtest.dir/flags.make
external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o: ../external/gmock/gtest/src/gtest-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/Rifdock/rifdock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o"
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gtest.dir/src/gtest-all.cc.o -c /home/longxing/Rifdock/rifdock/external/gmock/gtest/src/gtest-all.cc

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gtest.dir/src/gtest-all.cc.i"
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/Rifdock/rifdock/external/gmock/gtest/src/gtest-all.cc > CMakeFiles/gtest.dir/src/gtest-all.cc.i

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gtest.dir/src/gtest-all.cc.s"
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/Rifdock/rifdock/external/gmock/gtest/src/gtest-all.cc -o CMakeFiles/gtest.dir/src/gtest-all.cc.s

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires:

.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides: external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires
	$(MAKE) -f external/gmock/gtest/CMakeFiles/gtest.dir/build.make external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides.build
.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides

external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.provides.build: external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o


# Object files for target gtest
gtest_OBJECTS = \
"CMakeFiles/gtest.dir/src/gtest-all.cc.o"

# External object files for target gtest
gtest_EXTERNAL_OBJECTS =

external/gmock/gtest/libgtest.a: external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o
external/gmock/gtest/libgtest.a: external/gmock/gtest/CMakeFiles/gtest.dir/build.make
external/gmock/gtest/libgtest.a: external/gmock/gtest/CMakeFiles/gtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/Rifdock/rifdock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libgtest.a"
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && $(CMAKE_COMMAND) -P CMakeFiles/gtest.dir/cmake_clean_target.cmake
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/gmock/gtest/CMakeFiles/gtest.dir/build: external/gmock/gtest/libgtest.a

.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/build

external/gmock/gtest/CMakeFiles/gtest.dir/requires: external/gmock/gtest/CMakeFiles/gtest.dir/src/gtest-all.cc.o.requires

.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/requires

external/gmock/gtest/CMakeFiles/gtest.dir/clean:
	cd /home/longxing/Rifdock/rifdock/build/external/gmock/gtest && $(CMAKE_COMMAND) -P CMakeFiles/gtest.dir/cmake_clean.cmake
.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/clean

external/gmock/gtest/CMakeFiles/gtest.dir/depend:
	cd /home/longxing/Rifdock/rifdock/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/Rifdock/rifdock /home/longxing/Rifdock/rifdock/external/gmock/gtest /home/longxing/Rifdock/rifdock/build /home/longxing/Rifdock/rifdock/build/external/gmock/gtest /home/longxing/Rifdock/rifdock/build/external/gmock/gtest/CMakeFiles/gtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/gmock/gtest/CMakeFiles/gtest.dir/depend
