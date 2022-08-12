Mixing Many-Body and Pairwise-Additive Force Fields
===============
The purpose of this example is to demonstrate how to combine different molecules
using radically different pair pair styles together in the same simulation
(using moltemplate).


### Images

<img src="images/watMW.jpg" width=30> <img src="images/plus.svg" height=80> <img src="images/cyclopentane.jpg" width=80> <img src="images/rightarrow.svg" height=80>  <img src="images/cyclododecane+watMW_t=0ps_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80>  <img src="images/cyclododecane+watMW_t=50ps_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80>  <img src="images/cyclododecane+watMW_t=400ps_LR.jpg" width=150>

### Video

https://www.youtube.com/watch?v=IIIHg2p7QN4


### Details

This is a relatively complex example containing two different types of
coarse-grained molecules (water and cyclododecane)
which phase-separate during the simulation.
*(Note:These are united-atom models which do not have explicit hydrogen atoms.)*

1) This simulation uses the 3-body
(non-pairwise-additive, directional) coarse-grained "mW" water model:
Molinero, V. and Moore, E.B., J. Phys. Chem. B 2009, 113, 4008-4016
Simulations using the "mW" water model can be several orders of magnitude
faster than simulations using simple all-atom models such as SPCE or TIP3P.

2) A traditional united-atom TraPPE force field (which is pairwise-additive)
was used for the cyclododecane molecules.


Any force-field available in LAMMPS can be used with moltemplate.
New force-field styles are added by end users regularly.
For a current list, see:
https://docs.lammps.org/Commands_pair.html


### Instructions 
More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

### Step 1) README_setup.sh
This file explains how to use moltemplate.sh to build the files that
LAMMPS needs.

### Step 2) README_run.sh
This file explains how to use LAMMPS to run a simulation using the
files you created in step 1.


### Prerequisites

This example requires that LAMMPS was built with the "MANYBODY" package.
If lammps complains of a missing pair style, you will have to recompile
LAMMPS with the "MANYBODY" package enabled.
For details see:
https://lammps.sandia.gov/doc/Build_package.html
