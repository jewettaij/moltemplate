WARNINGS:

   This example was only intented to be a technical demonstration to show how
to combine totally different kinds of coarse-grained molecules (with different
kinds of force-fields) together in the same simulation in LAMMPS.  Tuning the
force-field parameters to get realistic results was not the goal.

1) If my understand is correct, this lipid model was not originally designed
to be used in vesicles with such high curvature.  It is not clear whether it
will behave realistically when placed in such a small vesicle.  The inner
layer shows signs of being in the gel phase at 300K, so for the example,
I raised the temperature to 345K to reduce this.  I suspect other lipid
models could behave differently at these temperatures and curvatures.

2) The protein model was originally globular and was heavily modified to
make it stable in a lipid bilayer.  Its behavior in a membrane is has
not been characterized.

   In addition, I have noticed that newer versions of PACKMOL do not
always succeed at generating a spherical vesicle in a reasonable amount of time.
(You may have to play with the .inp files in the packmol_files directory
 to get PACKMOL to produce any files at all.)
