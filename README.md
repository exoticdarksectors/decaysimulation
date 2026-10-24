# Summary
The folder decaysimulation contains the entirety of meson production and DY production pipelines, these include the meson generation methods via Pythia found in the mesongen folder and the decay of subsequent mesons into mCPs in the decay folder. Drell-Yan production is modelled in MadGraph, where the folder madgraph_mcp contains the Minimal_MCP model used for simulation. 

The mesongen program uses pythia -- a program for the generation of high-energy physics collision events -- to simulate the collision between two particles that can be set using the parameter files beam.config and momentum.config. mesonGen is an executable that simulates the collision and outputs a distribution of the momentum of the final state particles in the form of root files.

The decay program then takes these root files as inputs and saves the momentum distribution of mcp into another root file called mcp-production.root which represents a histogram of the following produced particles.

The madgraph_mcp folder contains a self-contained MadGraph setup (Minimal_MCP) that simulates the Drell–Yan production channel for millicharged particles (MCPs).


# Installation
Tested on macOS Ventura 13.6.7

C++ compiler info by performing `g++ --version`:

```
Apple clang version 14.0.3 (clang-1403.0.22.14.1)
Target: x86_64-apple-darwin22.6.0
Thread model: posix
```

## Dependencies

### ROOT and Cmake
The CERN's ROOT library, Pythia and the Cmake utility for building C/C++ code are necessary to compile the mesonGen and decayPion executables.

To do install ROOT and Cmake, execute the following commands in a terminal:

```
brew install root
brew install cmake
```
Define bash environment variables by adding the following lines to your `~/.bash_profile`:
```
export ROOTSYS=/usr/local/Cellar/root/6.32.06

export PATH=$PATH:$ROOTSYS/bin

```
for Linux users add:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
```
for MacOS users add:
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib
```

### Pythia
To install Pythia:

1. Go to https://www.pythia.org/
2. Download the latest version of Pythia as a tar ball, e.g. https://www.pythia.org/download/pythia83/pythia8312.tgz
3. Un-tar the archive:

```
tar -xvzf pythia8312.tgz
cd pythia8312
./configure
make -j8
```

Define bash environment variables by adding the following lines to your `~/.bash_profile`:

```
export PYTHIA8=/Users/leobailloeul/Documents/coding/software/pythia8309

export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc

export DYLD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PYTHIA8}/lib
```

where `PYTHIA8` needs to reflect the location where the Pythia library was compiled.

### MadGraph
To install Madgraph:

1. Go to http://madgraph.phys.ucl.ac.be/
2. Download the latest version of Madgraph as tar ball
3. Un-tar and place MadGraph somewhere in your `$HOME`:

```
# Example: install MG5_aMC v3.5.3 under ~/software
mkdir -p ~/software
tar -xzf MG5_aMC_v3.5.3.tar.gz
mv MG5_aMC_v3_5_3 ~/software/madgraph
```
Define bash environment variables by adding the following lines to your `~/.bash_profile`:

```
export MG5=$HOME/software/madgraph
export PATH=$PATH:$MG5/bin
```

### ExRootAnalysis (Used for MadGraph)

To install ExRootAnalysis:
1. Go to https://madgraph.mi.infn.it/Downloads/ExRootAnalysis/
2. Download the latest version of ExRootAnalysis as tar ball
3. Un-tar the archive:
```
tar -zxf ExRootAnalysis.tar.gz
cd ExRootAnalysis
make
```
Please refer to https://madgraph.mi.infn.it/Downloads/ExRootAnalysis/README in case of additional issues.

### CLion IDE
The CLion IDE is recommended for developing and debugging C/C++ code.
To install CLion, follow this [link](https://www.jetbrains.com/clion/download/#section=mac).

Request student/university free license by filling JetBrains' form.


# mesonGen

## Compilation
Edit the CMakeLists file and adjust the variable `IMPORTED_LOCATION` to the location where your Pythia library is installed: `IMPORTED_LOCATION "/Users/leobailloeul/Documents/coding/software/pythia8309/lib/libpythia8.dylib"`

### Clion IDE
1. Open CLion IDE
2. Select `Existing Project` mesonGen
3. Let CLion import the project by detecting CMakeLists.txt
4. To compile the project, select Build --> Buid Project


### Terminal
1. Open a bash terminal and go to the root of this github repo.
2. Create a build folder of your choice, e.g. `mkdir my-build`
3. Move to the build folder: `cd my-build`
4. Execute cmake: `cmake ..`
5. Build: `cmake --build .`
6. The executable is available in the `my-build` folder.


## Run
### Clion IDE
1. Edit the Run config, and pass parameters `seed` `core` `number of runs`
2. To run the decayPion executable, either to Run --> `Run mesonGen`, or click on the green play triangle on the upper right corner.

### Terminal

```
./mesonGen <seed> <core> <number of runs>
```

Example:

```
./mesonGen 1234 4 10000
```

Output files will be written int the `output-data` folder. One can combine output rootfiles using the command `hadd`.
TODO: make the output folder a parameter of the mesonGen binary.

# decayPion

## Compilation

### Clion IDE
1. Open CLion IDE
2. Select `Existing Project` decayPion
3. Let CLion import the project by detecting CMakeLists.txt
4. To compile the project, select Build --> Buid Project


### Terminal
1. Open a bash terminal and go to the root of this github repo.
2. Create a build folder of your choice, e.g. `mkdir my-build`
3. Move to the build folder: `cd my-build`
4. Execute cmake: `cmake ..`
5. Build: `cmake --build .`
6. The executable is available in the `my-build` folder.


## Run
### Clion IDE
1. Edit config, and pass parameters `full-path-to-input-root-file` `full-path-to-output-root-file`
2. To run the decayPion executable, either to Run --> `Run decayPion`, or click on the green play triangle on the upper right corner.

### Terminal

```
./decayPion <full-path-to-input-root-file> <full-path-to-output-root-file>
```

Example:

```
./decayPion /Users/leobailloeul/documents/coding/decaysimulation/mesongen/output-data/mesons_seed2024.root /Users/leobailloeul/documents/coding/decaysimulation/decay/output-data/mcp-production.root

```
# madgraph_mcp

1. Minimal_MCP: model for madgraph, minimal mcp lagrangian has been implemented in a way that can be feed into madgraph (search UFO, universial feynrules output? for more info). Move this into models folder in the madgraph folder
2. makeSetMchiMass: python script for making automation to loop through all mass points and generate madgraph events. Remeber to modify first lin into launch (your folder name) -n mass1
3. extractCross: read crossx.html and extract the cross section
4. getAgeo.cc: root macro that extracts the geoemetrical acceptance. x and y is the particles location at the detector. Modify dumplength.
