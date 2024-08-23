# CMakeLists.txt for the "event" package. It creates a library and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
project(decayPionProject)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT 6.28 CONFIG REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# Create a shared library.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.

# Create the main program using the library.
add_executable(decayPion decayPion.cc)
add_executable(decayVectorMeson decayVectorMeson.cc)
add_executable(2BodyMomentum decayVectorMeson.cc)

target_link_libraries(decayPion PUBLIC ROOT::Tree ROOT::Physics ROOT::TreePlayer)
target_link_libraries(decayVectorMeson PUBLIC ROOT::Tree ROOT::Physics ROOT::TreePlayer)
