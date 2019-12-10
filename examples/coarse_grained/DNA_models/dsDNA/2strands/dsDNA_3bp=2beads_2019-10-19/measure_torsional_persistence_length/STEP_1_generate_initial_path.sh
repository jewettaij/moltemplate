#!/usr/bin/env bash

DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
NMONOMERS=32
LMONOMER=1.0   #<-relaxed distance between monomers along the polymer axis (nm).
               #  (This is just a guess. We will have to run some simulations at
               #   the desired temperature to learn what the actual length is.)

for ((i=0; i<$NMONOMERS; i++)); do 
    XCOORD=0.0
    YCOORD=0.0
    ZCOORD=`echo "$LMONOMER*($i+0.5-0.5*$NMONOMERS)" | bc -l`
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

