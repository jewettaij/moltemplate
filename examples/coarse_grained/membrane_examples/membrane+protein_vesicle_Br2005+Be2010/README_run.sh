# --- Running LAMMPS ---
#  -- Prerequisites: --
# The "run.in.nvt" file is a LAMMPS input script containing
# references to the input scripts and data files
# you hopefully have created earlier with MOLTEMPLATE and PACKMOL:
#   system.in.init, system.in.settings, system.in.coords, system.data,
#   and table_int.dat
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# Assuming "lmp_mpi" is the name of the LAMMPS binary,
# run lammps in this order:

lmp_serial -i run.in.min  # Minimize the system (important and very slow)

lmp_serial -i run.in.make_uniform  # Trap the lipids between concentric
                             # spherical shells and equilibrate.  This insures
                             # that the lipids are distributed uniformly
                             # on the spherical surface.  (Unfortunately,
                             # PACKMOL does not guarantee this.)

lmp_serial -i run_T=345K.in # Let the vesicle relax and equilibrate.
                            # (Temperature was chosen to avoid the membrane
                            #  gel phase which occurs with this model in
                            #  vesicles of this size.)



#If you have compiled the MPI version of lammps, you can run lammps in parallel
# mpirun -np 16 lmp_mpi -i run.in.min
#   or
# mpirun -np 16 lmp_mpi -i run.in.nvt
# (assuming you have 16 processors available)
