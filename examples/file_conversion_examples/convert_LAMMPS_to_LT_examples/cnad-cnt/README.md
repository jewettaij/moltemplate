Combining, editing, and duplicating LAMMPS DATA files
==========================================
This example uses the *ltemplify.py* program (included with moltemplate) to convert LAMMPS data files (and input scripts) created by other 3rd-party programs into MOLTEMPLATE-compatible .lt files.  These files can be used to construct a new system using the components of the original system.  In step 1, *ltemplify.py* is used to extract two different molecules (CNAT, CNT) the [cnat-cnt.data](cnat-cnt.data) and [cnat-cnt.in](cnat-cnt.in) files.  In step 2, *moltemplate.sh* is used to make copies of these molecules and to move and rotate them.  In the process a new input script and DATA file will be generated.

The *ltemplify.py* program is explained [here](../../../../doc/doc_ltemplify.md).

### Images
Original system:

<img src="images/cnad-cnt_orig.jpg" width=180>

Molecules extracted by *ltemplify.py*:

<img src="images/cnad.jpg" height=100> <img src="images/plus.svg" height=80> <img src="images/cnt.jpg" height=100>

New system created by *moltemplate.sh*:

<img src="images/cnad-cnt_after_rotate_copy.jpg" width=200>


#### Disclaimer
The molecules in this example are (intentionally) not physically realistic.  The purpose of this example is to demonstrate ltemplify usage.


### Instructions
Run the commands (follow the instructions) in these files:

1) [README_step1_run_ltemplify.sh](README_step1_run_ltemplify.sh).  This will create two files: *cnad.lt*, *cnt.lt* which define the *"CNAD"* and *"CNT"* molecules.  These files are referenced in [system.lt](system.lt).  

2) [README_step2_run_moltemplate.sh](README_step2_run_moltemplate.sh).  This creates new LAMMPS data and input files: system.data, system.in.init, system.in.settings.  (These files are referenced in [run.in.nvt](run.in.nvt), used in step 3.)

3) OPTIONAL: Run a LAMMPS simulation

To run a short LAMMPS simulation, you can use the [run.in.nvt](run.in.nvt) file, for example:
```
lmp_mpi -i run.in.nvt
```
#### Note
The name of the LAMMPS binary ("lmp_mpi") may vary

#### Note
**Because the original force field parameters were intentionally altered, the system will move in an unrealistic way.**  *The goal of this example is only to demonstrate how to use "ltemplify.py" to convert LAMMPS input and data files into LT format and back again.*


### Details

#### Required input files

1) [cnad-cnt.data](cnad-cnt.data)
This is a LAMMPS data file containing the coordinates and the topology for a system combining the two molecules together.  ltemplify will extract molecules from this file, one at a time.

2) [cnad-cnt.in](cnad-cnt.in)
 This file contains force-field parameters and old run settings for the system.  (We ignore the run settings in this file.)  The force-field parameters in the "cnad-cnt.in" file are only necessary because we are going to build a completely new set of simulation input files. (We are not only going to rotate them and duplicate the molecules.)  ltemplify.py will extract the force field parameters from this file.  This approach allows us to combine these molecules with other types of molecules later on.)

3) [system.lt](system.lt)
The "system.lt" contains the instructions what we will do with these molecules after ltemplify.py has converted them into .LT format.  In this example it contains instructions for rotating and copying the two molecules, (It also defines the periodic boundary conditions.)


#### Credit
The original version of this example was contributed by Aysun Itai.
