#!/usr/bin/env bash

# The following lines of BASH script generate a 3-column coordinate file
# ($DEST_FILE) containing the coordinates for a circular ring of points.
# Later we will use "genpoly_lt.py" to generate a polymer by placing
# a monomer at every one of these points.
# This can be done more easily using moltemplate's .move() and .rot() commands.
# But in this example, I wanted to show how to create a polymer that
# follows the shape of a general curve.

DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
N_MONOMERS=100
L_MONOMER=13.761  #<--relaxed distance between monomers along the polymer axis
PI=3.141592653589793
R_CIRCLE=`echo "$L_MONOMER*$N_MONOMERS/(2*$PI)" | bc -l`
for ((i=0; i<$N_MONOMERS; i++)); do 
    XCOORD=`echo "$R_CIRCLE*c(($i+0.5)*2*$PI/$N_MONOMERS)" | bc -l`
    YCOORD=`echo "$R_CIRCLE*s(($i+0.5)*2*$PI/$N_MONOMERS)" | bc -l`
    ZCOORD=0.0
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

