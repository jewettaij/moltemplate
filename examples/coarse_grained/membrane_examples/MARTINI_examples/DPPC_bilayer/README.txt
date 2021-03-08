
This is an example of a coarse-grained DPPC lipid-bilayer (MARTINI model).
Lipid orientation within the membrane was chosen randomly.
(This was done using the "new random" command to select from versions of the
 lipid molecule which were rotated by different amounts around the Z axis.
 See the "system.lt" file for details.)
It's probably a good idea to run the simulation for a few ns to allow
the the orientation and position of the lipids to relax.

step 1)
To build the files which LAMMPS needs, follow the instructions in:
README_setup.sh

step 2)
To run LAMMPS with these files, follow these instructions:
README_run.sh

------- CITE -----------------------------
NOTE: We extracted the parameters in the MARTINI force field from the files
distributed with the "EMC" tool.  If you use these .lt files, please also cite:
P. J. in ‘t Veld and G. C. Rutledge, Macromolecules 2003, 36, 7358.
