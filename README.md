# Gravity-Simulation

This repository contains scripts that allow you to run a physical simulation of a system with gravity.
The same principles could however be used to simulate any forces you like.
The time integration is done using the Verlet algorithm.

The Jupyter Notebook contains explanations and motivations for this implementation.
On top of that there is a Python script that can run the simulation in one go for convenience.
The same simulation has also been implemented in C++, for which three files are relevant: 
the source code script itself, the header file that defines parameters and initial conditions, 
and a Python script to plot the results.
Any outputs from the scripts are placed in the 'output' directory.
The Python script outputs two plots, for the trajectories in the xy-plane and the energies in the system.
The C++ script outputs a csv file with the results, which is used as input for the Python plotting script, 
which then creates the same two plots as the Python simulation script.

For the test that was run and can be seen in the output directory here, 
the results from both the Python and C++ simulation are the same, as they should.
The only difference is the runtime; on my machine the Python simulation took 3.93 seconds, 
while the C++ simulation only needed 0.96 seconds.
