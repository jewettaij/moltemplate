#!/usr/bin/env bash


DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"
NPOINTS=32
#LMONOMER=1.0108   #<--relaxed distance between monomers along the polymer axis
LMONOMER=1.00       #<--I shortened it a little to start with a smaller circle
                    #   (Be sure to run minimization beforehand to relax 
                    #    the bonds to their rest lengths.)
PI=3.141592653589793
RCIRCLE=`echo "$LMONOMER*$NPOINTS/(2*$PI)" | bc -l`
for ((i=0; i<$NPOINTS; i++)); do 
    XCOORD=`echo "$RCIRCLE*c(($i+0.5)*2*$PI/$NPOINTS)" | bc -l`
    YCOORD=`echo "$RCIRCLE*s(($i+0.5)*2*$PI/$NPOINTS)" | bc -l`
    ZCOORD=0.0
    echo "$XCOORD $YCOORD $ZCOORD" >> "$DEST_FILE"
done

