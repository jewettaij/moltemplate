# PREREQUISITES:
#    system.data, system.in.init, system.in.settings
# To create these files, follow the instructions in README_setup.sh
#
# Then, run LAMMPS using the following command:

lmp_mpi -i run.in

# NOTE: To run this in parallel, use:
# mpiexec -np 4 lmp_mpi -i run.in

# NOTE: The name of the LAMMPS binary typically varies depending upon how it
# was compiled.  If the name of your LAMMPS binary is different
# (eg "lmp_ubuntu"), replace "lmp_mpi" with the correct binary name.

