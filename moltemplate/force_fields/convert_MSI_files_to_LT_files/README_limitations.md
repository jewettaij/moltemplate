The "msifrc2lt.py" program was originally written to convert FRC files
(storing force field parameters in ACCELRYS FRC format)
into moltemplate LT files.  This conversion program was used to convert
the "compass_published.frc" file into "compass_published.lt".

Unfortunately, "msifrc2lt.py" currently only supports the COMPASS force field.
The PCFF and CVFF force fields are much more complicated and do not work
with "msifrc2lt.py".  There are several reasons for this (explained at the
end of this message).

Most people who use these force fields are using ACCELRYS software
(such as Materials Studio).  But this software is very expensive.

Unfortunately, as of 2022-3-19, I have no plans to fix "msifrc2lt.py".


# Recommendations:

If you need to use either PCFF or CVFF, try using EMC
https://LAMMPS.sandia.gov/prepost.html
EMC supports the PCFF, CVFF, and COMPASS force fields.
Can determine the correct atom types and looking up force field parameters
for these force fields.  (I don't know if it understands "auto-equivalences".)
The main limitation of EMC is that you don't have very much control over the
position of the molecules in your simulation.  EMC will try to place the
molecules in realistic locations randomly.  If this is not what you want,
you can use EMC to create a LAMMPS DATA file and INPUT SCRIPT for a
simulation containing a SINGLE COPY of the molecule(s) you want to simulate.
Then you can use moltemplate to make many copies of this molecule
(or molecules) and arrange them in space manually where you want each
molecule to go.

1) Download EMC from http://montecarlo.sourceforge.net
2) Unpack it.
3) Optional: If you have an alternate version of PCFF file
(such ass "pcff_iff.frc"), then copy your custom file to the
"field/pcff/" subdirectory and rename the new file "pcff.frc"
4) Run EMC.

The "ltemplify.py" program is explained here:
https://github.com/jewettaij/moltemplate/blob/master/doc/doc_ltemplify.md

Other resources to consider:
(I have never tried these resources, so I don't know if they work.)

a) msi2lmp
https://github.com/lammps/lammps/tree/master/tools/msi2lmp
(WARNING: The original author of this code is no longer maintaining it, 
and people have complained in the mailing list that it does not work
very reliably.)

b) struct2lammps
https://nanohub.org/resources/struct2lammps/


If you are curious, I decided not to fix the problems with "msifrc2lt.py"
because:

- Both the CVFF and PCFF force fields use "auto equivalences" to lookup
bonded interactions between atoms whenever the normal atom force-field
lookup rules fail..  "auto equivalences" are very complicated.  Last I checked,
the "msi2lmp" tool completely ignores "auto equivalences", which I think is a
serious flaw.  (This means that your resulting data file might be missing
angles, dihedrals, or impropers which should be present in your molecule.
Strangely, most people don't seem bothered by this.)  Moltemplate does not
ignore "auto equivalences", but it turned out to be too difficult to
implement these correctly given the time I had, so I gave up.  Fixing
this limitation is difficult.  Others have tried it and given up.
It might be easier to start from scratch than it would be to try
editing the "msifrc2lt.py" file.

- The COMPASS, PCFF, and CVFF force fields are proprietary.
The only files I have found for these force fields are either incomplete
(eg "compass_published.frc") or encrypted.  This means even if I were to
fix the limitations in "msifrc2lt.py" to make the conversion possible,
only people who have access to an unencrypted version of the 
"pcff.frc" or "cvff.frc" files would be able to take advantage of this feature.
Alternatively, the incomplete versions of these force fields
(such as "compass_published.frc") lack atom types (such as SP2 carbons)
needed to build most of the molecules that people care about.
For me, it wasn't worth the considerable effort needed to get "msifrc2lt.py"
to work with these force fields.

-Andrew 2022-3-20
