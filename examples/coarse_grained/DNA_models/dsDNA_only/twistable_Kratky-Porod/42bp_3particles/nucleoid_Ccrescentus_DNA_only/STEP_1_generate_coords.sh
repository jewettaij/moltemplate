# Create a file containing the coordinates for the locations of each monomer
# of a circular polymer stretched out into a straight conformation
# with N_MONOMERS_HALF monomers going in each direction.
# Later we will place monomers (in our polymer model) at each of these locations
#
#          *--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
#          |                                                     |
#          *--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
#
# Because I am doing this in BASH instead of a proper programming language
# I decided it was easier to do this in two steps:
# First use integer coordinates.
# Later scale them by the size of each monomer.
# Feel free to implement this some other way.
# First we have to define how many monomers there are and the spacing between
# them.

N_MONOMERS_HALF=47619 # The number of monomers in each direction.  (The total
                      # number monomers in the circle is twice this number.)

L_MONOMER=13.944     # physical distance between monomers along the polymer axis

rm -f crds.raw
for (( i=0; i < N_MONOMERS_HALF; i++ )); do
    echo "$i 0.5 0.0" >> crds.raw
done
for (( i=N_MONOMERS_HALF-1; i>= 0; i-- )); do
    echo "$i -0.5 0.0" >> crds.raw
done

# Now scale the coordinates by the size of each monomer (L_MONOMER), and center
# them (by subtracting half of the length of the ring from the x coordinate),
# and save the resulting file in "init_crds_polymer_backbone.raw".

awk -v L=$L_MONOMER -v N=$N_MONOMERS_HALF \
    '{printf("%.9g %.8g %.8g\n", ($1-0.5*N)*L, $2*L, $3*L)}' < crds.raw \
    > moltemplate_files/init_crds_polymer_backbone.raw

# Clean up:
rm -f crds.raw
