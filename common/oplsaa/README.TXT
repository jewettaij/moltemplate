Credit: OPLSAA force-field provided by Jason Lambert.
        (Thanks! -Andrew)

This directory contains a sample moltemplate file ("oplsaa.lt") containing
force-field definitions relevant to the system you are simulating.
It should contain force-field information for a subset of the atoms in the
oplsaa forcefield.  This directory also contains the python code 
("oplsaa_moltemplate.py") used to generate that file, as well as the full
oplsaa force field file: "oplsaa.txt".

python oplsaa_moltemplate.py  <force field file name>

First, back-up the oplsaa.txt file and then delete all the unnecessary atoms
from the file.

It is ideal for you to run this on a subset of the OPLS forcefield relevant 
to your problem.  (It is possible for you to generate a full OPLS force field
moltemplate file but that demands a lot of time and generates a file that is
on the order 100 gigabytes, which moltemplate can not read.)  
Save yourself time and energy, and make a copy of the oplsaa.txt that only
contains atoms relevant to your problem.

Let's assume this new file is called "oplsaa_subset.txt".  If you were 
working with ethylene you would delete every atom except for atom 85(C=) 
and 86(H-C=) and use that as the input file for oplsaa_moltemplate.py
You do not need to make any other changes other than that.  Then run:

python oplsaa_moltemplate.py oplsaa_subset.txt

  This will create a file named "oplsaa.lt"

Look over the newly created "oplsaa.lt" file.  Then copy it to the "common" directory (where other force-field files like gaff.lt are located, ...or just copy it to wherever you plan to run moltemplate).  

Then you can create molecules that refer to it using:

Ethylene inherits OPLSAA {
  # Definition of Ethylene molecule goes here.
}

CHARGE:

Currently, you must assign the atom charges yourself, but changes to
moltemplate may be added so that this is not necessary.  For now, use the set
command in LAMMPS:

http://lammps.sandia.gov/doc/set.html
For example:
set type 3 charge 0.5  # Note this number (3) may not match the OPLS atom number
                       # Check the output_ttree/ttree_assignments.txt file.

-Jason Lambert
February, 2014

Please email bugs to jewett.aij@gmail.com

