Random Heteropolymer Example
====================

This is an example of how to build a polymer out of randomly-chosen monomers.
In this case, monomers will be chosen at random from two types
(denoted "2bead" and "3bead", although you can have as many types as you like).
You can also constrain the end-caps to be a particular type (eg "3bead").

#### Images

<img src="images/2bead.jpg" height=50> <img src="images/plus.svg" height=80>
<img src="images/3bead.jpg" height=40> <img src="images/rightarrow.svg" height=80>
<img src="images/random_heteropolymer_30_20_t=0.jpg" width=300> <img src="images/rightarrow.svg" height=80>
<img src="images/random_heteropolymer_30_20_t=700ps.jpg" width=180> 

Note that you still must manually connect the monomers together with bonds.
(by adding linse to the "Data Bond List" section of the "polymer.lt" file
located in the moltemplate_files directory).

Also note that the properties of the bonds connecting monomers
(ie length, rigidity) will be automatically determined,
depending on atoms at that location in the polymer.
Similarly 3-body angle, and 4-body dihedral interactions
will be determined automatically
(according to rules inside the "forcefield.lt" file).

For more details, see the "Random Arrays" section of the moltemplate manual.
*([here](https://moltemplate.org/doc/moltemplate_manual.pdf#subsection.8.4) and
[here](https://moltemplate.org/doc/moltemplate_manual.pdf#subsubsection.8.4.1))*


### Instructions 
Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

#### Step 1) README_setup.sh
This file explains how to use moltemplate.sh to build the files that
LAMMPS needs.

#### Step 2) README_run.sh
This file explains how to use LAMMPS to run a simulation using the
files you created in step 1.
