This directory contains the analysis scripts and instructions that would be
needed to measure the torsional persistence length of the DNA polymer

First copy the ../moltemplate_files directory into this directory:

cp -r ../moltemplate_files .

Then run the following scripts:

./STEP_1_generate_initial_path.sh
./STEP_2_generate_LAMMPS_files.sh
./STEP_3_run_sim.sh     # <- this file needs editing to match your environment
./STEP_4_measure_torsional_persistence_length.sh

The torsional persistence length should be printed to the terminal,
as well as other measurements (length-per-bp, twist-per-bp, etc...).

