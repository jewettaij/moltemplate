# Here we generate the starting coordinates of the simulation
# using PACKMOL.



# You must run each packmol commend one after the other.
#
# NOTE: If PACKMOL gets stuck in an endless loop, then edit the corresponding
# "inp" file (edit the "tolerance" and "distscale").  This should not happen,
# but keep in mind that PACKMOL can be very slow.  (Step3 takes several days.)
#
# NOTE: You can also usually interrupt packmol after 1 hour.  The solution
#       at that point (an .xyz file) should be good enough for use.
#      (Of course, you should check to see that the file exists and the
#       structure looks reasonable visually before interrupting packmol.
#       Small holes are okay. The simulation protocol we use should close them.)

packmol < step1_proteins.inp   # This step determines the protein's location
                               # It takes ~40 minutes (on a 2015 intel i7)
packmol < step2_innerlayer.inp # This step builds the inner monolayer
                               # It takes ~10 hours
packmol < step3_outerlayer.inp # This step builds the outer monolayer
                               # It takes ~1-3 days


# NOTE: PLEASE USE "packmol", NOT "ppackmol".  ("ppackmol" is the
#       parallel-version of packmol using OpemMP.  This example has NOT been
#       tested with "ppackmol".  Our impression was that the "ppackmol"
#       version is more likely to get stuck in an infinite loop. -Andrew 2015-8)
