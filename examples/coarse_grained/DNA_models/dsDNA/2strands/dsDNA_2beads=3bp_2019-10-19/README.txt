# This example demonstrates how to prepare a simulation of double-stranded (ds)
# DNA using the coarse-grained "3bp=2beads" model.
#
#    FEATURES:
# 1)This model has a double-helical conformation with a major and minor grooves.
#    It can be used to model interaction of DNA with simple protein models
#    that bind specifically to the major or minor grooves.
# 2) The polymer can be bent arbitrarily far without kinking or numerical
#    instability.  This allows the use of comparatively large timesteps.
# 3) The masses were adjusted to make it possible to set the timestep to "1.0".
#
#    (NOTE: Use of a common timestep makes it easier to combine different coarse
#     grained models together in the same simulation.  Using heavier masses can
#     compensate for stiff bonds or bond-angles.  Using unrealistic masses is
#     justified if the processes-of-interest occur on time-scales far larger
#     than the vibrational spectra of individual bonds or angles in the coarse
#     grained model.  At these time scales, the rates of diffusive motion can
#     be controlled by adjusting the Langevin parameters and solvent properties,
#     and are not effected by the mass of the particles of the model.)
#
# Geometrical and mechanical properties of this DNA model (simulated at T=300K)
# are compared with experiment below:
#                                                |  this model  | experimenal |
# ----------------------------------------------------------------------------
# persistence length:                            |  48.7101 nm  |     50 nm   |
# torsional persistence length:                  | 109.9752 nm  |    110 nm   |
# helical twist angle between monomers (3bp):    | 102.7797 deg |  see below  |
# -> helical twist angle between base-pairs:     |  34.2599 deg |   34.3 deg  |
# distance along the axis between monomers (3bp):|  0.98293 nm  |  see below  |
# -> distance along the axis between base-pairs: |  0.32764 nm  |   0.332 nm  |
# ----------------------------------------------------------------------------
#
#    PREREQUISITES:
# LAMMPS must be compiled with the "MOLECULE" AND "USER-MISC" packages enabled.
# If you receive this error message (or something similar):
#    dihedral_style spherical: Unknown dihedral style
# you must recompile LAMMPS. (https://lammps.sandia.gov/doc/Build_package.html)
#
#    WARNING:
# These files (originally uploaded on 2019-12-10) contain many comments
# (beginning with "#") which are probably misleading and no longer relevant.
# I will try to clean up these files over time.
# Please let me know if anything doesn't work.
#
#    INSTRUCTIONS:
#
# Instructions on how to build LAMMPS input files and 
# run a short simulation are provided below
#
# The following file contain instructions explaining how to generate
# the curve that you want the DNA polymer to follow.
# (You can also run it as an executable.)

./STEP_1_generate_initial_path.sh

# The next file explains how to convert this curve into a moltemplate file, and
# how to run moltemplate on that file. (You can also run it as an executable.)

./STEP_2_generate_LAMMPS_files.sh

# Finally, to run the LAMMPS simulation follow the instructions in this file:
# STEP_3_run_sim.sh
# You will have to edit the file to specify the name of the LAMMPS binary
# you are using (for example, "lmp_ubuntu"), and the number of processors
