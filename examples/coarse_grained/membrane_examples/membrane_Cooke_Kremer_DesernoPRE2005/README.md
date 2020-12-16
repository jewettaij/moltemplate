Cooke & Kremer (2005) coarse grained membrane
=============================================
This example contains an implementation of the DPPC lipid bilayer described in: "Tunable generic model for fluid bilayer membranes" Cooke IR, Kremer K, Deserno M, Physical Review E, 2005.  *(That paper describes phase separation of lipids in multicomponent membranes.  This simulation includes only one type of lipid.  However there is another moltemplate example featuring multicomponent membranes.)

### Images

<img src="images/CDlipid.jpg" height=100> <img src="images/rightarrow.svg" height=80> <img src="images/CDlipid_bilayer_t=0_nopbc_occ.jpg" width=190>  <img src="images/rightarrow.svg" height=80> <img src="images/CDlipid_bilayer_t=600000steps_npt_occ.jpg" width=180>


### Strategy

First, a hexagonal array of lipids is created.  Rectangular periodic boundary conditions are applied, and the system is equilibrated at zero tension and constant temperature.  To change the number of molecules in the system and initial simulation box size, edit the [system.lt](moltemplate_files/system.lt) file.  To change the simulation conditions, edit the [run.in.npt](run.in.npt) file.


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)
