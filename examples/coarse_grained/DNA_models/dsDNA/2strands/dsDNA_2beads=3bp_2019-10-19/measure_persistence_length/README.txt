This directory contains the analysis scripts and instructions that would be
needed to measure the persistence length of the DNA polymer
(as well as other geometrical properties you might be interested in).

First copy the ../moltemplate_files and ../run* scripts into this directory:

cp -r ../moltemplate_files .
cp ../run*.in .

Then run the following scripts:

./STEP_1_generate_initial_path.sh
./STEP_2_generate_LAMMPS_files.sh
./STEP_3_run_sim.sh     # <- this file needs editing to match your environment
./STEP_4_measure_distances_and_angles.sh    # (optional)
./STEP_5_measure_persistence_length.sh

The persistence length should be printed to the terminal
(...in monomers, not in nm.  To convert to nm, measure the length-per-monomer by
 running the scripts in the "measure_torsional_persistence_length" directory.)
