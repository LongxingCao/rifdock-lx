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
CMAKE_SOURCE_DIR = /home/longxing/devel/rifdock-lx

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/longxing/devel/rifdock-lx/build

# Include any dependencies generated for this target.
include external/gmock/CMakeFiles/gmock_main.dir/depend.make

# Include the progress variables for this target.
include external/gmock/CMakeFiles/gmock_main.dir/progress.make

# Include the compile flags for this target's objects.
include external/gmock/CMakeFiles/gmock_main.dir/flags.make

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o: external/gmock/CMakeFiles/gmock_main.dir/flags.make
external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o: ../external/gmock/gtest/src/gtest-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifdock-lx/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o -c /home/longxing/devel/rifdock-lx/external/gmock/gtest/src/gtest-all.cc

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifdock-lx/external/gmock/gtest/src/gtest-all.cc > CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifdock-lx/external/gmock/gtest/src/gtest-all.cc -o CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.requires:

.PHONY : external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.requires

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.provides: external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.requires
	$(MAKE) -f external/gmock/CMakeFiles/gmock_main.dir/build.make external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.provides.build
.PHONY : external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.provides

external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.provides.build: external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o


external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o: external/gmock/CMakeFiles/gmock_main.dir/flags.make
external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o: ../external/gmock/src/gmock-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifdock-lx/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/src/gmock-all.cc.o -c /home/longxing/devel/rifdock-lx/external/gmock/src/gmock-all.cc

external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock-all.cc.i"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifdock-lx/external/gmock/src/gmock-all.cc > CMakeFiles/gmock_main.dir/src/gmock-all.cc.i

external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock-all.cc.s"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifdock-lx/external/gmock/src/gmock-all.cc -o CMakeFiles/gmock_main.dir/src/gmock-all.cc.s

external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.requires:

.PHONY : external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.requires

external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.provides: external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.requires
	$(MAKE) -f external/gmock/CMakeFiles/gmock_main.dir/build.make external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.provides.build
.PHONY : external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.provides

external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.provides.build: external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o


external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: external/gmock/CMakeFiles/gmock_main.dir/flags.make
external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: ../external/gmock/src/gmock_main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifdock-lx/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.o -c /home/longxing/devel/rifdock-lx/external/gmock/src/gmock_main.cc

external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock_main.cc.i"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifdock-lx/external/gmock/src/gmock_main.cc > CMakeFiles/gmock_main.dir/src/gmock_main.cc.i

external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock_main.cc.s"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifdock-lx/external/gmock/src/gmock_main.cc -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.s

external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.requires:

.PHONY : external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.requires

external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.provides: external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.requires
	$(MAKE) -f external/gmock/CMakeFiles/gmock_main.dir/build.make external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.provides.build
.PHONY : external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.provides

external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.provides.build: external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o


# Object files for target gmock_main
gmock_main_OBJECTS = \
"CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o" \
"CMakeFiles/gmock_main.dir/src/gmock-all.cc.o" \
"CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"

# External object files for target gmock_main
gmock_main_EXTERNAL_OBJECTS =

external/gmock/libgmock_main.a: external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o
external/gmock/libgmock_main.a: external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o
external/gmock/libgmock_main.a: external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o
external/gmock/libgmock_main.a: external/gmock/CMakeFiles/gmock_main.dir/build.make
external/gmock/libgmock_main.a: external/gmock/CMakeFiles/gmock_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/devel/rifdock-lx/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libgmock_main.a"
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean_target.cmake
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/gmock/CMakeFiles/gmock_main.dir/build: external/gmock/libgmock_main.a

.PHONY : external/gmock/CMakeFiles/gmock_main.dir/build

external/gmock/CMakeFiles/gmock_main.dir/requires: external/gmock/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o.requires
external/gmock/CMakeFiles/gmock_main.dir/requires: external/gmock/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o.requires
external/gmock/CMakeFiles/gmock_main.dir/requires: external/gmock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.requires

.PHONY : external/gmock/CMakeFiles/gmock_main.dir/requires

external/gmock/CMakeFiles/gmock_main.dir/clean:
	cd /home/longxing/devel/rifdock-lx/build/external/gmock && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean.cmake
.PHONY : external/gmock/CMakeFiles/gmock_main.dir/clean

external/gmock/CMakeFiles/gmock_main.dir/depend:
	cd /home/longxing/devel/rifdock-lx/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/devel/rifdock-lx /home/longxing/devel/rifdock-lx/external/gmock /home/longxing/devel/rifdock-lx/build /home/longxing/devel/rifdock-lx/build/external/gmock /home/longxing/devel/rifdock-lx/build/external/gmock/CMakeFiles/gmock_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/gmock/CMakeFiles/gmock_main.dir/depend

