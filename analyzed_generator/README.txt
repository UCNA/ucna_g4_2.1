March 15, 2016
Xuan Sun

The code in this directory is designed to be run after GEANT4 finishes a simulation.
It turns the output .txt file into a .root TTree which has more useful branches.
This was designed originally to remove the ROOT libraries from the GEANT4 CMake compilation process.
Now ROOT compiles successfully, but in the past I had a lot of problems with compiling both.

Also now it isn't perfectly sensible since ROOT is incorporated in PrimaryGeneratorAction. But whatever.

Moral of the story: GEANT4 gives you a UCNASimOutput.txt (or whatever you choose to call it).
You compile ttree.cc and ttree.hh using 'make' and execute it.
In ttree.cc, you'll provide the name (possibly path too) of the simulation output file.
