Coarse grained membrane (Brannigan et al. 2005)
==================
This example contains an implementation of the DPPC lipid bilayer described in:
1) G. Brannigan, P.F. Philips, and F.L.H. Brown, Physical Review E, Vol 72, 011915 (2005)
2) M.C. Watson, E.S. Penev, P.M. Welch, and F.L.H. Brown, J. Chem. Phys. 135, 244701 (2011)
A truncated version of this lipid (named "DLPC") has also been added.

### Images

<img src="images/DPPC.jpg" height=100> <img src="images/plus.svg" height=80> <img src="images/DLPC.jpg" height=82> <img src="images/rightarrow.svg" height=80> <img src="images/DPPC+DLPC_bilayer32x37_t=0ps_no_pbc_LR.jpg" width=190>  <img src="images/rightarrow.svg" height=80> <img src="images/DPPC+DLPC_bilayer32x37_t=500ps_LR.jpg" width=180>


### Strategy
First, a hexagonal array of lipids is created.  Rectangular periodic boundary conditions are applied, and the system is equilibrated at zero tension and constant temperature.  To change the composition of the membrane, edit the [system.lt](moltemplate_files/system.lt) file.  To change the simulation conditions, edit the [run.in.npt](run.in.npt) file.


### Prerequisites

1) This example requires custom features to be added to LAMMPS, which usually require LAMMPS to be compiled from source code.  If you encounter the error *"Invalid pair_style"*, then you must add these features to LAMMPS.  To do that:
a) Download the LAMMPS source code (if you have not yet done so), either from https://lammps.sandia.gov/download.html, or using *"git clone https://github.com/lammps/lammps ~/lammps_src"*
b) download the "additional_lammps_code" from http://moltemplate.org (upper-left corner menu)
c) unpack it
d) copy the .cpp and .h files to the src folding of your lammps installation.
e) (re)compile LAMMPS.

2) This example uses "dihedral_style fourier" which requires a version of LAMMPS compiled with support for the optional "EXTRA-MOLECULE" package.  If you encounter the error *"Invalid dihedral_style"*, then see [this page](https://lammps.sandia.gov/doc/Build_package.html) for instructions to compile LAMMPS to support this package.


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)


### Details

The DPPC lipid bilayer described in:
1) G. Brannigan, P.F. Philips, and F.L.H. Brown, Physical Review E, Vol 72, 011915 (2005)
2) M.C. Watson, E.S. Penev, P.M. Welch, and F.L.H. Brown, J. Chem. Phys. 135, 244701 (2011)

As in Watson(JCP 2011), rigid bond-length constraints, have been replaced by harmonic bonds.

A truncated version of this lipid (named "DLPC") has also been added.  The bending stiffness of each lipid has been increased to compensate for the additional disorder resulting from mixing two different types of lipids together.  (Otherwise pores appear.) Unlike the original "DPPC" molecule model, the new "DPPC" and "DLPC" models have not been carefully parameterized to reproduce the correct behavior in a lipid bilayer mixture.


#### Known issues
This is not an ideal coarse grained lipid mixture.  Simulations of small liposomes with high curvature have shown that the the DPPC model at 300K is likely crystalize and phase separate from the DLPC model lipids.
