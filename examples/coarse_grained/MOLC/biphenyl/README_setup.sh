# The moltemplate files we need are in the "moltemplate_files" directory.

# Enter that directory

cd moltemplate_files/

  # Then run moltemplate.sh to generate the files LAMMPS needs.
  #
  # NOTE: This example uses an unconventional particle representation
  # ("hybrid molecular ellipsoid").  Consequently moltemplate should 
  # be run with the "-atomstyle" argument to inform it where to find
  # the attributes of each atom.  (Incidentally, this list of attributes 
  # can be overridden and replaced with an arbitrary list attributes.)
  # For this model, the "-overlay bonds" argument is also needed, allowing
  # multiple bonds between the same pair of atoms.  (Moltemplate allows this.)

 moltemplate.sh -atomstyle "hybrid molecular ellipsoid" -overlay-bonds system.lt

  # This will create the following files which LAMMPS needs:
  #    system.data system.in.settings system.in.init

  # Copy them back to the parent directory where we will run LAMMPS:

  mv -f system.data system.in.settings system.in.init ../

  # Optional: Cleanup the files we no longer need

  rm -rf output_ttree system.in

cd ../

