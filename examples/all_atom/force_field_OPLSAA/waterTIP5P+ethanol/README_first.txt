This example of a mixture of rigid TIP5P water
with an ordinary flexible molecule (ethanol)
was contributed by Cheng Chen (Queen Mary University) on 2018-8-02.

The "fix rigid" command (or "fix rigid/nph/small") was used to integrate 
the equations of motion for the TIP5P molecule.
...and a combination of "fix nve" was used to integrate
the equations of motion for the flexible ethanol molecule.
In this example "fix langevin" was used to regulate temperature.

However there are multiple other ways to run these kinds of simulations.
A combination of "fix rigid/npt/small" and "fix nvt" might also work
(without "fix langevin").  For a discussion of alternative ways to simulate
mixtures of flexible and rigid molecules, see the "run.in.npt" file
(and the "WARNING.txt" file).

Follow the instructions in "README_setup.sh" to create the files which LAMMPS
needs to start the simulation.  The instructions for running LAMMPS are in
"README_run.sh"

(The instructions in "README_remove_irrelevant_info.sh" are optional.)
