# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /u/group/clas12/packages/cmake/3.15.2/bin/cmake

# The command to remove a file.
RM = /u/group/clas12/packages/cmake/3.15.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build

# Include any dependencies generated for this target.
include CMakeFiles/LowEnergy.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LowEnergy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LowEnergy.dir/flags.make

CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o: CMakeFiles/LowEnergy.dir/flags.make
CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o: ../LowEnergyReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o -c /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/LowEnergyReader.cpp

CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/LowEnergyReader.cpp > CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.i

CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/LowEnergyReader.cpp -o CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.s

# Object files for target LowEnergy
LowEnergy_OBJECTS = \
"CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o"

# External object files for target LowEnergy
LowEnergy_EXTERNAL_OBJECTS =

LowEnergy: CMakeFiles/LowEnergy.dir/LowEnergyReader.cpp.o
LowEnergy: CMakeFiles/LowEnergy.dir/build.make
LowEnergy: libclas12reader.a
LowEnergy: libhipochain.a
LowEnergy: /u/apps/root/6.18.04/root/lib/libPhysics.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libPostscript.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libROOTDataFrame.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libROOTVecOps.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libRint.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libTreePlayer.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libTree.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libEG.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libGraf3d.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libGpad.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libGraf.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libHist.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libMatrix.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libMathCore.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libMultiProc.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libImt.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libNet.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libRIO.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libThread.so.6.18.04
LowEnergy: /u/apps/root/6.18.04/root/lib/libCore.so.6.18.04
LowEnergy: CMakeFiles/LowEnergy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LowEnergy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LowEnergy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LowEnergy.dir/build: LowEnergy

.PHONY : CMakeFiles/LowEnergy.dir/build

CMakeFiles/LowEnergy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LowEnergy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LowEnergy.dir/clean

CMakeFiles/LowEnergy.dir/depend:
	cd /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build /work/clas12/users/esteejus/clas12analysisForFarm/MagneticFieldStudy/build/CMakeFiles/LowEnergy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LowEnergy.dir/depend

