Ethylene, Benzene mixture using PACKMOL
==============
This example also shows how to use moltemplate in together with PACKMOL to create mixture of two small organic molecules using the *OPLSAA* force field.  [PACKMOL](http://m3g.iqm.unicamp.br/packmol/home.shtml) was used to generate the initial coordinates of the system.  PACKMOL is a useful program for generating atomic coordinates.  In this example, moltemplate.sh is only used to create the topology, force-field and charges, and PACKMOL generates the coordinates, which moltemplate reads (in "step 1").  Moltemplate can also be used for generating atomic coordinates, especially for mixing many small molecules together.  However I wanted to demonstrate how to combine PACKMOL with moltemplate.sh.  In some other scenarios, such as protein solvation, PACKMOL does a much better job than moltemplate.


#### Images
<img src="images/ethylene.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/benzene.jpg" width=110> <img src="images/plus.svg" height=80> PACKMOL <img src="images/rightarrow.svg" height=80> <img src="images/ethylene+benzene_box80x80x80_LR.jpg" width=200>

The simulation size and number of ethylene and benzene is specified in the [mix_ethylene+benzene.inp](./packmol_files/mix_ethylene+benzene.inp) and [system.lt](./moltemplate_files/system.lt) files.  (The numbers in these files must agree.)  The simulation contitions can be controlled by editing the [run.in.npt](run.in.npt) file.

## Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)

