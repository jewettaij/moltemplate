#!/usr/bin/env bash


DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
NMONOMERS=32
LMONOMER=1.0   #<-relaxed distance between monomers along the polymer axis (nm).
               #  (This is just a guess. We will have to run some simulations at
               #   the desired temperature to learn what the actual length is.)
PI=3.141592653589793
RCIRCLE=`echo "$LMONOMER*$NMONOMERS/(2*$PI)" | bc -l`
for ((i=0; i<$NMONOMERS; i++)); do 
    XCOORD=`echo "$RCIRCLE*c(($i+0.5)*2*$PI/$NMONOMERS)" | bc -l`
    YCOORD=`echo "$RCIRCLE*s(($i+0.5)*2*$PI/$NMONOMERS)" | bc -l`
    ZCOORD=0.0
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

