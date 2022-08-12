polymer translocation
=====================

This example contains a (crude and somewhat simple) example of
the translocation of a (rather short) polymer through a hole in a wall,
surrounded by an explicit LJ solvent.

### Images

<img src="images/polymer_LR.jpg" width=170> <img src="images/plus.svg" height=80> <img src="images/walls_LR.jpg" width=150> <img src="images/plus.svg" height=80> <img src="images/solvent_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80> <img src="images/walls+solvent+polymer_t=0.jpg" height=200>

#### Video

https://www.youtube.com/watch?v=hxAufspeg2s


### Instructions 
Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

#### Step 1) README_setup.sh
This file explains how to use moltemplate.sh to build the files that
LAMMPS needs.

#### Step 2) README_run.sh
This file explains how to use LAMMPS to run a simulation using the
files you created in step 1.

Note that there are two ways to run this simulation, at constant volume
or at constant pressure.
When the simulation is run at constant pressure,
the polymer and the liquid is expelled through
the openning as the simulation box expands.
*(Running simulations containing immobile objects at constant
pressure in LAMMPS is complicated.  See the "run.in.npt" file for details.)*


### Prerequisites

To run the simulation at constant pressure
LAMMPS must be compiled with the "RIGID" package enabled.
So if LAMMPS generates the following error:
"rigid: Unknown fix", then you must follow
[these instructions](https://lammps.sandia.gov/doc/Build_package.html),
and recompile LAMMPS.
