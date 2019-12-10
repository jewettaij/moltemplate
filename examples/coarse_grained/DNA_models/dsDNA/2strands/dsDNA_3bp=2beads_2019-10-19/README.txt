# This example demonstrates how to build a simulation of coarse-grained
# double-stranded ("ds") DNA using the "3bp=2beads" model.
# This model has a double-helical conformation with a major and minor groove.
# The polymer can be bent arbitrarily far without causing numeric instability.
# It's geometrical and mechanical properties are compared with experiment below:
#                                             | simulated    | experimenal |
# --------------------------------------------------------------------------
# persistence length:                         |  48.7101 nm  |     50 nm   |
# torsional persistence length:               | 109.9752 nm  |    110 nm   |
# helical twist angle between base-pairs:     |  34.2599 deg |   34.3 deg  |
# distance along the axis between base-pairs: |  0.32764 nm  |  0.332 nm   |
# --------------------------------------------------------------------------
#
# Instructions on how to build LAMMPS input files and 
# run a short simulation are provided below

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
