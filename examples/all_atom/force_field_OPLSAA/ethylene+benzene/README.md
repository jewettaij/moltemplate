Ethylene, Benzene mixture
==============
A mixture of two small organic molecules using the *OPLSAA* force field.  In this example, the ethylene molecules were initially arranged in a rectangular lattice.  The benzene molecules were also arranged in a lattice, and were shifted to avoid overlap with the ethylene molecules. *(See the [system.lt file](./moltemplate_files/system.lt) for details.)*

<img src="images/ethylene.jpg" width=110>
<img src="images/plus.svg" height=80>
<img src="images/benzene.jpg" width=110>
<img src="images/rightarrow.svg" height=80>
<img src="images/ethylene+benzene_t=0_LR.jpg" width=150>
<img src="images/rightarrow.svg" height=80>
<img src="images/ethylene+benzene_50bar_t=100000_LR.jpg" width=150>

*(Alternatively, you can create a single lattice and specify the number of ethelene and benzene molecules you want in it using moltemplate's "new random([],[])" command, which is explained in the manual.  You can also use PACKMOL to create random mixtures of molecules.)*

## Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*