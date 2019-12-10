#!/usr/bin/env bash


DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
NMONOMERS=32
LMONOMER=0.98293   #<--relaxed distance between monomers along the polymer axis
PI=3.141592653589793
RCIRCLE=`echo "$LMONOMER*$NMONOMERS/(2*$PI)" | bc -l`
for ((i=0; i<$NMONOMERS; i++)); do 
    XCOORD=`echo "$RCIRCLE*c(($i+0.5)*2*$PI/$NMONOMERS)" | bc -l`
    YCOORD=`echo "$RCIRCLE*s(($i+0.5)*2*$PI/$NMONOMERS)" | bc -l`
    ZCOORD=0.0
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

