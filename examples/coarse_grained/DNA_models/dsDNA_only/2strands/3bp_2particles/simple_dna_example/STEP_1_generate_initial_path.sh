#!/usr/bin/env bash


DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
N_MONOMERS=250
L_MONOMER=0.98293   #<--relaxed distance between monomers along the polymer axis
PI=3.141592653589793
R_CIRCLE=`echo "$L_MONOMER*$N_MONOMERS/(2*$PI)" | bc -l`
for ((i=0; i<$N_MONOMERS; i++)); do 
    XCOORD=`echo "$R_CIRCLE*c(($i+0.5)*2*$PI/$N_MONOMERS)" | bc -l`
    YCOORD=`echo "$R_CIRCLE*s(($i+0.5)*2*$PI/$N_MONOMERS)" | bc -l`
    ZCOORD=0.0
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

