#!/usr/bin/env bash

DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
NPOINTS=32
LMONOMER=1.0108   #<--relaxed distance between monomers along the polymer axis

for ((i=0; i<$NPOINTS; i++)); do 
    XCOORD=0.0
    YCOORD=0.0
    ZCOORD=`echo "$LMONOMER*($i+0.5-0.5*$NPOINTS)" | bc -l`
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

