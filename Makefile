# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /mnt/c/Users/Davide/Desktop/AG-FIT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Davide/Desktop/AG-FIT

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/c/Users/Davide/Desktop/AG-FIT/CMakeFiles /mnt/c/Users/Davide/Desktop/AG-FIT/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/c/Users/Davide/Desktop/AG-FIT/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named FitGA

# Build rule for target.
FitGA: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FitGA
.PHONY : FitGA

# fast build rule for target.
FitGA/fast:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/build
.PHONY : FitGA/fast

src/Chromosome.o: src/Chromosome.cc.o

.PHONY : src/Chromosome.o

# target to build an object file
src/Chromosome.cc.o:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Chromosome.cc.o
.PHONY : src/Chromosome.cc.o

src/Chromosome.i: src/Chromosome.cc.i

.PHONY : src/Chromosome.i

# target to preprocess a source file
src/Chromosome.cc.i:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Chromosome.cc.i
.PHONY : src/Chromosome.cc.i

src/Chromosome.s: src/Chromosome.cc.s

.PHONY : src/Chromosome.s

# target to generate assembly for a file
src/Chromosome.cc.s:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Chromosome.cc.s
.PHONY : src/Chromosome.cc.s

src/DataGenerator.o: src/DataGenerator.cc.o

.PHONY : src/DataGenerator.o

# target to build an object file
src/DataGenerator.cc.o:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/DataGenerator.cc.o
.PHONY : src/DataGenerator.cc.o

src/DataGenerator.i: src/DataGenerator.cc.i

.PHONY : src/DataGenerator.i

# target to preprocess a source file
src/DataGenerator.cc.i:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/DataGenerator.cc.i
.PHONY : src/DataGenerator.cc.i

src/DataGenerator.s: src/DataGenerator.cc.s

.PHONY : src/DataGenerator.s

# target to generate assembly for a file
src/DataGenerator.cc.s:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/DataGenerator.cc.s
.PHONY : src/DataGenerator.cc.s

src/FitGA.o: src/FitGA.cc.o

.PHONY : src/FitGA.o

# target to build an object file
src/FitGA.cc.o:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/FitGA.cc.o
.PHONY : src/FitGA.cc.o

src/FitGA.i: src/FitGA.cc.i

.PHONY : src/FitGA.i

# target to preprocess a source file
src/FitGA.cc.i:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/FitGA.cc.i
.PHONY : src/FitGA.cc.i

src/FitGA.s: src/FitGA.cc.s

.PHONY : src/FitGA.s

# target to generate assembly for a file
src/FitGA.cc.s:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/FitGA.cc.s
.PHONY : src/FitGA.cc.s

src/ParametersDomain.o: src/ParametersDomain.cc.o

.PHONY : src/ParametersDomain.o

# target to build an object file
src/ParametersDomain.cc.o:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/ParametersDomain.cc.o
.PHONY : src/ParametersDomain.cc.o

src/ParametersDomain.i: src/ParametersDomain.cc.i

.PHONY : src/ParametersDomain.i

# target to preprocess a source file
src/ParametersDomain.cc.i:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/ParametersDomain.cc.i
.PHONY : src/ParametersDomain.cc.i

src/ParametersDomain.s: src/ParametersDomain.cc.s

.PHONY : src/ParametersDomain.s

# target to generate assembly for a file
src/ParametersDomain.cc.s:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/ParametersDomain.cc.s
.PHONY : src/ParametersDomain.cc.s

src/Population.o: src/Population.cc.o

.PHONY : src/Population.o

# target to build an object file
src/Population.cc.o:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Population.cc.o
.PHONY : src/Population.cc.o

src/Population.i: src/Population.cc.i

.PHONY : src/Population.i

# target to preprocess a source file
src/Population.cc.i:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Population.cc.i
.PHONY : src/Population.cc.i

src/Population.s: src/Population.cc.s

.PHONY : src/Population.s

# target to generate assembly for a file
src/Population.cc.s:
	$(MAKE) -f CMakeFiles/FitGA.dir/build.make CMakeFiles/FitGA.dir/src/Population.cc.s
.PHONY : src/Population.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... FitGA"
	@echo "... edit_cache"
	@echo "... src/Chromosome.o"
	@echo "... src/Chromosome.i"
	@echo "... src/Chromosome.s"
	@echo "... src/DataGenerator.o"
	@echo "... src/DataGenerator.i"
	@echo "... src/DataGenerator.s"
	@echo "... src/FitGA.o"
	@echo "... src/FitGA.i"
	@echo "... src/FitGA.s"
	@echo "... src/ParametersDomain.o"
	@echo "... src/ParametersDomain.i"
	@echo "... src/ParametersDomain.s"
	@echo "... src/Population.o"
	@echo "... src/Population.i"
	@echo "... src/Population.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

