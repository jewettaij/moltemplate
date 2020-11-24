Macroscopic example
===========================
This is a simulation of the collapse of a stack of spherical particles arranged into a pyramid shape.  An initial (small) disturbance in the shape of the pyramids leads to their eventaul collapse in an avalanche.

### Images

<img src="images/pyramids_vs_gravity_t=04800steps_LR.jpg" width=300>
<img src="images/rightarrow.svg" height=80>
<img src="images/pyramids_vs_gravity_t=12200steps_LR.jpg" width=300>
<img src="images/rightarrow.svg" height=80>
<img src="images/pyramids_vs_gravity_t=33000steps_LR.jpg" width=300>


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files:

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)


### Details

The particles experience a downward force similar to gravity.  They roll down the pyramids and bounce off the "ground". The bouncing is due to a repulsive external force which is placed there using LAMMPS "fix wall/lj126" feature.  (See the [run.in](run.in) file.)  To aid in visualization, an immobile hexagonal array of particles is also placed on the bottom surface as scenery.  (Those particles does not exert any force on the moving particles.)

In this example, the cannon balls are initially slightly too far apart from each other.  (The system is not in a local energy minimum.)  As the simulation proceeds, they will begin to move slightly.  The resulting movement of the balls propogates through the pyramids, eventually causing an avalanche approximately 5000 timesteps later.  (If I had minimized the system first, or added frictional damping, the pyramids would not collapse.)

LAMMPS can be used to run [realistic simulations of granular flow](https://lammps.sandia.gov/doc/Howto_granular.html), however this simulation does not take advantage of this capability.  This simulation contains only simple point-like Lennard-Jones particles.  There is no friction in this simulation.  *(The particle heights should eventually approach the Boltzmann distribution for some temperature which is consistent with the initial gravitational energy of the system.)*
