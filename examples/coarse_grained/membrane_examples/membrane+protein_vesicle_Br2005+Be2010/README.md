Lipid membrane vesicle with membrane protein inclusions
====================
This example shows how to build a spherical vesicle made from DPPC lipids and 120 trans-membrane protein inclusions.  The coordinates for the ingredients of the vesicle are generated by PACKMOL (a 3rd-party tool).  This is a slow example because PACKMOL has difficulty packing molecules uniformly on a surface.  Several simulations are necessary to close holes in the membrane created by PACKMOL.

### Images

<img src="images/DPPC.jpg" height=100> <img src="images/plus.svg" height=80> <img src="images/4HelixBundle_LR.jpg" height=130> <img src="images/rightarrow.svg" height=80>
<img src="images/vesicle_membrane+protein_LR.jpg" width=400>


### Prerequisites

1) This example requires PACKMOL.  You can download PACKMOL [here](http://www.ime.unicamp.br/~martinez/packmol/)  (Moltemplate does not come with an easy way to generate spherically-symmetric structures, so I used the PACKMOL program to move the lipids into position.)

2) This example requires custom features to be added to LAMMPS, which usually require LAMMPS to be compiled from source code.  If you encounter the error *"Invalid pair_style"*, then you must add these features to LAMMPS.  To do that:
a) Download the LAMMPS source code (if you have not yet done so), either from https://lammps.sandia.gov/download.html, or using *"git clone https://github.com/lammps/lammps ~/lammps_src"*
b) download the "additional_lammps_code" from http://moltemplate.org (upper-left corner menu)
c) unpack it
d) copy the .cpp and .h files to the src folding of your lammps installation.
e) (re)compile LAMMPS.

3) This example uses "dihedral_style fourier", which requires a version of LAMMPS compiled with support for the optional "USER-MISC" package.  If you encounter the error *"Invalid dihedral_style"*, then see [this page](https://lammps.sandia.gov/doc/Build_package.html) for instructions to compile LAMMPS with the "USER-MISC" package enabled.


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)


### Details

This example contains a coarse-grained model of a 4-helix bundle protein inserted into a lipid bilayer (made from a mixture of DPPC and DLPC).


#### Protein Model

The coarse-grained protein is described in:
G. Bellesia, AI Jewett, and J-E Shea, Protein Science, Vol19 141-154 (2010)
Here we use the "AUF2" model described in that paper.  (In the paper, the protein is soluble and the hydrophobic beads face inwards.  In this simulation, the protein is embedded in a lipid bilayer and the hydrophobic beads face outwards towards the lipids.)


#### Memebrane Model

The DPPC lipid bilayer described in:
1) G. Brannigan, P.F. Philips, and F.L.H. Brown, Physical Review E, Vol 72, 011915 (2005)
2) M.C. Watson, E.S. Penev, P.M. Welch, and F.L.H. Brown, J. Chem. Phys. 135, 244701 (2011)

As in Watson(JCP 2011), rigid bond-length constraints, have been replaced by harmonic bonds.


--- Building the files necessary to run a simulation in LAMMPS ---

step 1) Run PACKMOL

        Type these commands into the shell.
        (Each command could take several hours.)

cd packmol_files
  packmol < step1_proteins.inp   # This step determines the protein's location
  packmol < step2_innerlayer.inp # this step builds the inner monolayer
  packmol < step3_outerlayer.inp # this step builds the outer monolayer
cd ..

step 2) Run MOLTEMPLATE
        Type these commands into the shell.
        (This could take up to 10 minutes.)

cd moltemplate_files
  moltemplate.sh system.lt -xyz ../system.xyz
  mv -f system.in* system.data ../
  cp -f table_int.dat ../
cd ..

--- Running LAMMPS ---

step3) Run LAMMPS:
        Type these commands into the shell.
        (This could take days.)

lmp_mpi -i run.in.min  # Minimize the system (important, and very slow)
lmp_mpi -i run.in.make_uniform
lmp_mpi -i run_T=360K.in  # Run a simulation at constant tmperature

If you have compiled the MPI version of lammps, you can run lammps in parallel:

mpiexec -np 4 lmp_mpi -i run.in.min
mpiexec -np 4 lmp_mpi -i run.in.make_uniform
mpiexec -np 4 lmp_mpi -i run.in_T=360K

(Assuming you have 4 cores, for example.)