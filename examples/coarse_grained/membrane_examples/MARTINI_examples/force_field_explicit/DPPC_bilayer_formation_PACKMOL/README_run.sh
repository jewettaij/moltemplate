# --- Running LAMMPS ---
# -------- PREREQUISITES: --------
# The 2 files "run.in.min", "run.in.anneal", and "run.in.nvt" are LAMMPS
# input scripts which link to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_mpi -i run.in.min     # minimization
lmp_mpi -i run.in.anneal  # high->low temp annealing simulation to form the
                          # bilayer.  By the end of the simulation, the
                          # system is now at T=300K, pressure=1bar

# There is no guarantee that a lipid bilayer has formed.
# Be sure to check that you have a smooth closed lipid bilayer before
# proceeding.  (To check what the membrane looks like, follow the instructions
# in the "README_visualization" file.)
# If the lipids are not in the shape of a well-formed membrane,
# then change the PACKMOL random seed, rebuild the system, and run the
# simulation again.  (For details, see the "README.txt" file.)

lmp_mpi -i run.in.nvt     # simulation at constant volume (optional)

# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.npt
#mpirun -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
