You can visualize the trajectory (dump file) with the program OVITO
(www.ovito.org). You will have to define few parameters, though.
Invoke ovito using:

ovito traj.lammpstrj

On the right panel of ovito, under the section "LAMMPS dump"
click on "File contains a time series" to display the entire trajectory.

Then click the button "Edit column mapping" and enter the fields as specified
in the png file located here:

https://github.com/jewettaij/moltemplate/blob/master/examples/coarse_grained/ellipsoids_CG_benzene/README_visualization_OVITO_ellipsoids.png

(However, in this particular example, the "Orientation" columns are named
 "c_quat" instead of "c_q".  As before, they are ordered: W, X, Y, Z.)

Note: The polymer will not appear as it does in the oxDNA litterature.
      Instead each nucleic acid will appear as a single ellipsoid.
Why:
      In the LAMMPS implementation of oxDNA, quaternions are used to keep track
      of the orientation of each nucleic acid so that LAMMPS can internally
      determine the positions of multiple internal particles located within
      each nucleic acid.  If my understanding is correct, then these particles
      are held rigid relative to each other during the simulation.  In order
      to associate quaternions with each nucleic acid, the LAMMPS data file
      which describes the position of the nucleic acids in the polymer
      represents each nucleic acid as a single ellipsoidal shaped particle.
      (That's why it uses "atom_style hybrid bond ellipsoid").  Internally
      however, LAMMPS will behave as if there are multiple particles per
      nucleic acid and calculate the forces and torques on it accordingly.
      OVITO knows nothing about these internal particles, and will simply
      show the nucleic acid as an ellipsoid whose orientation reflects
      the orientation of that nucleic acid.
