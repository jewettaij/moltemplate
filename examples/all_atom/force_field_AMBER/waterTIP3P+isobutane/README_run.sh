# --- Running LAMMPS ---
# -------- REQUIREMENTS: ---------
# This example requires building LAMMPS with the "EXTRA-MOLECULE" package.
# (because it makes use of "gaff.lt" which uses dihedral_style fourier).
# If you encounter the error *"Invalid dihedral_style"*, then this page:
# https://docs.lammps.org/Build_package.html
# for instructions to compile LAMMPS to support this package.
# -------- PREREQUISITES: --------
# The 2 files "run.in.npt", and "run.in.nvt" are LAMMPS
# input scripts which link to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_mpi -i run.in.npt  # minimization and simulation at constant pressure
lmp_mpi -i run.in.nvt  # simulation at constant volume

#(Note: The constant volume simulation lacks pressure equilibration. These are
#       completely separate simulations. The results of the constant pressure
#       simulation might be ignored when beginning the simulation at constant
#       volume.  (This is because restart files in LAMMPS don't always work,
#       and I was spending a lot of time trying to convince people it was a
#       LAMMPS bug, instead of a moltemplate bug, so I disabled restart files.)
#       Read the "run.in.nvt" file to find out how to use the "read_restart"
#       command to load the results of the pressure-equilibration simulation,
#       before beginning a constant-volume run.





# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.npt
#mpirun -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
