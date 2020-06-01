#!/usr/bin/env bash


DEST_FILE="moltemplate_files/init_crds_polymer_backbone.raw"
rm -f "$DEST_FILE"

N_MONOMERS=59105
N_IHF=1070
L_MONOMER=0.98293   # relaxed distance between monomers along the polymer
WIDTH_IHF=10        # how many monomers is IHF bound to?


# ---- Optional: Shorten the physical length of the polymer. ----
#           We do this because the polymer has N_IHF bend modifications, and
#           each one reduces the polymer length by roughly 2*WIDTH_IHF monomers.
#           Reducing the length of L_MONOMER will cause the monomers to be
#           closer together.  (We will minimize the system before the main
#           simulation to prevent steric collisions).  I think reducing the
#           overall polymer length will improve numerical stability, by not
#           overstretching the polymer at the beginning of the simulation.
#           (But I don't think it is strictly necessary to do this.)

L_MONOMER=`awk -v N=$N_MONOMERS -v L=$L_MONOMER -v n=$N_IHF -v w=$WIDTH_IHF 'BEGIN{print L*(N-n*w*2)/N}'`

# ---- (end of) Optional ----


# Now use a for-loop to create coordinates in a straight line

rm -f crds.raw
for (( i=0; i < N_MONOMERS; i++ )); do
    echo "$i 0 0" >> crds.raw
done

# This will create a file "crds.raw" containing:
# 0 0 0
# 1 0 0
# 2 0 0
# 3 0 0
# :

# Now center and re-scale these coordinates using awk
#    (Note: I pref awk to using BASH bc for floating point math.)
awk -v L=$L_MONOMER -v N=$N_MONOMERS \
    '{printf("%.9g %.8g %.8g\n", ($1-0.5*N)*L, $2*L, $3*L)}' < crds.raw \
    > moltemplate_files/init_crds_polymer_backbone.raw

rm -f crds.raw
