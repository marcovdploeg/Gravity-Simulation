# Gravity Simulation

This repository contains scripts that allow you to run a physical simulation of a system with gravity.
The same principles could however be used to simulate any forces you like.
The simulation also includes (inelastic) collisions, different boundary conditions and external forces.
The time integration is done using the (velocity) Verlet algorithm.

## Theory notebooks

The theory behind the simulation can be found in the "Theory" directory.
Here: 
- The Jupyter notebook "Gravity\_simulation\_theory.ipynb" contains explanations and motivations for
the basic version of this implementation. 
- The notebook "Collision\_simulation\_theory.ipynb" contains the explanations for 
the implementation of the collision physics.
- The markdown file "Derivations.md" contains some theoretical physical derivations for these collisions.
- The notebook "Boundary\_conditions\_theory.ipynb" contains the explanations for 
the implementation of the boundary conditions, consisting of the 'no/infinite box', 
'finite box' and 'periodic box' scenarios.
- Finally, the notebook "External\_force\_theory.ipynb" contains the explanations for the implementation
of external forces, that can be defined however and with whichever parameters one desires.
- An additional notebook "Many\_particles\_theory.ipynb" explains how to set up a simulation with many particles 
where the initial positions are uniformly distributed and the initial velocities are drawn from the 
Maxwell-Boltzmann distribution. Some improvements to the algorithm are also made here.

## Simulation scripts and output

On top of this, there is a Python script that can run the simulation in one go for convenience.
The same simulation has also been implemented in C++, for which three files are relevant: 
the source code script itself, the header file that defines parameters and initial conditions, 
and a Python script to plot the results.
Any outputs from the scripts are placed in the "example_output" directory.
The Python script outputs two plots, for the trajectories in the xy-plane and the energies in the system.
Additionally it creates a text file with the parameters that were used, for later reference.
The C++ script outputs a csv file with the results, which is used as input for the Python plotting script, 
which then creates the same two plots as the Python simulation script.

## Examples

For each theory notebook, corresponding test runs are included and can be seen in the "example_output" directory. 
One is a basic test, the second tests collisions and the third tests boundary conditions. 
The basic test shows some nice elliptical orbits of four objects around their centre of mass.
The test with collisions shows how four objects bounce around while losing energy, 
so that they end up clustering together.
The test with boundary conditions shows how objects will move differently under the periodic boundary conditions, 
which also apply the minimal image convention to calculate forces between objects.
The test with external forces shows how objects move under the influence of this external force while 
it dominstes over the gravitational force between the objects. In this case sinusoidal movement occurs.
The test with many particles shows the simulation is also capable of simulating e.g. gasses, 
although the runtime will grow large with more objects.

The results from the Python and C++ simulations are the same, as they should.
The only difference is the runtime; on my machine the Python simulation took about 1.2 seconds, 
while the C++ simulation only needed about 0.5 seconds for the basic example 
(and only 0.05 seconds when compiled for maximum speed).
