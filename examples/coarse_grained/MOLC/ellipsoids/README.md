## Summary

This is a simple example to demonstrate how to setup simulations containing
ellipsoids using moltemplate.  The ellipsoids interact with each other 
via the Gay-Berne potential.

## Motivation:

This is a preliminary example that was created before the paper describing
the MOLC coarse graining method was published.  The ellipsoids in this
example were originally intended to represent a coarse-grained version of
benzene molecules.  But do not expect realistic behaviour from this example.
The long-range electrostatic interactions are missing.

## Instructions:

To build the files which LAMMPS needs, follow the instructions in
README_setup.sh

To run the simulation in LAMMPS, follow the instructions in:
README_run.sh

Finally, to view the simulation results in OVITO, follow the instructions in:
- README_visualization_OVITO.txt
- README_visualization_OVITO_ellipsoids.png (in the parent directory)

## Credits:

This example was provided by Otello M. Roscion (U.Southampton)
and Matteo Ricci(U.Bologna)  Many thanks to them for editing and debugging
moltemplate code to get this working!
