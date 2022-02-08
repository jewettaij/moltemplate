Water, sodium, chloride mixture
==============
A mixture SPC/E water with sodium and chloride ions.  The water molecules were initially arranged in a rectangular lattice.  The sodium and chloride ions were also arranged in a lattice, and were shifted to avoid overlap with the water molecules and each other.  *(Alternatively, you can create a single lattice and specify the number of water, sodium, and chloride ions you want in it using moltemplate's "new random([],[])" command, which is explained in [the manual](https://moltemplate.org/doc/moltemplate_manual.pdf#subsubsection.8.9.1).  This gives you more control over the concentration of each ingredient.  You can also use PACKMOL to create random mixtures of molecules.)*


#### Images

<img src="images/wat.jpg" width=80> <img src="images/plus.svg" height=80> <img src="images/Na.jpg" width=80> <img src="images/plus.svg" height=80> <img src="images/Cl.jpg" width=80> <img src="images/rightarrow.svg" height=80>  <img src="images/waterSPCE+Na+Cl_t=0.jpg" width=150> <img src="images/rightarrow.svg" height=80>  <img src="images/waterSPCE+Na+Cl_t=100ps.jpg" width=150>

The number of water molecules and ions in the simulation and and the simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### Details 

Here I am using the SPC/E water model with long range electrostatics.  Monovalent ion parameters for Ewald and SPC/E water are from [Joung S, Cheatham TE, JPCB, 2008, 112(30):9020-41](https://doi.org/10.1021/jp8001614).  See the end of the [ions.lt](moltemplate_files/ions.lt) file for details.


### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*
