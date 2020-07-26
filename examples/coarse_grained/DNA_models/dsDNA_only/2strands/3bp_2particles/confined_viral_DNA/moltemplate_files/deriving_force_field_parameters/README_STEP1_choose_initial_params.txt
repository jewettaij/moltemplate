To generate a reduced-atom (a.k.a. "coarse-grained") model of DNA, we chose
a subset of the atoms in real DNA where we would place coarse-grained beads.

In particular, for this model, we chose to place coarse-grained beads
at the locations of the C3' atoms in the backbones of real DNA,
at EVERY THIRD BASE PAIR.

Then we added bonds (and also angle and dihedral interactions) to these beads
and chose their parameters (eg. spring stiffnesses) so that the resulting
simulated polymer had the same structure (and behavior) of real DNA.

To do this, the force-field parameters for the coarse grained DNA model were
chosen to reproduce the geometry and shape heterogeneity found in structures
of DNA found in the PDB database (at https://www.rcsb.org).
It was not feasible (or advisable) to measure these distances and angles
for every structure in the PDB database.  Instead, a (hopefully) unbiased
subset of PDB files containing DNA structures were selected using the
procedure explained here:

https://github.com/jewettaij/dlpdb/blob/master/examples/dna_example/README_STEP_1_prepare_pdb_files.sh

  --- UGLY DETAILS: ---
Then these files were processed to discarding overhangs, end-caps, and 
molecules which were not DNA.  (This process was semi-automated.  The
resulting PDB files were inspected visually and (occasionally) discarded.)
(Then, the order of the bases in the PDB file was changed ("interleaved") so
that the double-stranded DNA is treated as a single polymer (instead of 2
polymers), and the residue-ID/SeqNum numbers increase along the length of
that (double stranded) polymer.)
  ---------------------

Then the distrubution of distances and angles (between various combinations
of atoms from different bases) was calculated using the procedure outlined here:

https://github.com/jewettaij/dlpdb/blob/master/examples/dna_example/statistics_keeping_every_3rd_base_pair/calc_distances_angles.sh

These distributions were fit to the normal distribution (by calculating
the average and standard deviations of these measured quantities).

An initial guess for the force-field parameters was made in order to
reproduce the observed distribution of measured quantities.
For example, it was assumed that the distribution between a given pair of
bonded beades (eg. "r"), is given by a Boltzmann distribution:
   P(r) = A * exp(-0.5*k*(r-r0)^2 / (kB*T))
Matching these with the measured distribution of lengths from the PDB files
yeilds: r0 = average distance between that pair of atoms in the pdb file,
and k = kB*T / σ^2, where "σ" is the observed standard deviation, and "T" is
the temperature of the simulation we plan to run.
This was repeated for every bond type, (and angle type, and dihedral type)
which existed in the coarse grained model.

This crude approach provides a reasonable first guess for the force-field
parameters, but it ignores the effects of constraints and other nearby
atoms on these measured distances and angles.  So we will have to refine
these parameters later to compensate for these other effects.

There is a picture of the resulting coarse-grained DNA model here:

https://github.com/jewettaij/dlpdb
