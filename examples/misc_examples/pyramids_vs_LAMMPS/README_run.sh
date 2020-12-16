# --- Running LAMMPS ---
#
#  -- Prerequisites: --
# The "run.in" input script depends on the following files which
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_mpi -i run.in  # simulate the effects of gravity on your system


# Note: The name of your lammps binary (eg "lmp_mpi") may vary.

# If you have compiled the MPI version of lammps, you can run lammps in parallel
# mpirun -np 4 lmp_mpi -i run.in
# (assuming you have 4 processors available)
