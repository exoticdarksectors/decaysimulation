# CMakeLists.txt for the "event" package. It creates a library and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
project(mesonGenProject)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT 6.28 CONFIG REQUIRED)

# Find Pythia8 -- hardcoded paths for now
set(PYTHIA8 "/Users/leobailloeul/Documents/coding/software/pythia8309")
set(PYTHIA8_INCLUDE_DIR ${PYTHIA8}/include)
set(PYTHIA8_LIBRARY ${PYTHIA8}/lib)
set(PYTHIA8DATA ${PYTHIA8}/share/Pythia8/xmldoc)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${PYTHIA8_INCLUDE_DIR})

# Create a shared library.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.

# Create the main program using the library.
add_executable(mesonGen mesonGen.cc)

add_library(pythia8 SHARED IMPORTED) # or STATIC instead of SHARED
set_target_properties(pythia8 PROPERTIES
        IMPORTED_LOCATION "/Users/leobailloeul/Documents/coding/software/pythia8309/lib/libpythia8.dylib"
        INTERFACE_INCLUDE_DIRECTORIES ${PYTHIA8_INCLUDE_DIR}
)
target_link_libraries(mesonGen PUBLIC ROOT::Tree ROOT::Physics pythia8)