#!/usr/bin/env bash

#dump2data.py -raw -tstart 1 < traj.lammpstrj | ./raw2blockaverage.py 2 \
#             > traj_center-of-mass_trace.raw
#
# "-tstart 1" skips the first frame.  
# The simulation might need more than one dump interval of time before 
# equilibrating.  In that case, use a larger -tstart value.

TSTART=10000000  # How many verlet iterations of simulation do we throw away?

# Old method: Find the center-of-mass position of the atoms in each base, and
#             assume that the polymer axis lies along the displacement vector
#             connecting one center-of-mass position with the next one.
#             THIS DOES NOT WORK, because the center-of-mass is not very
#             close to the center of the polymer axis.
#dump2data.py -raw -tstart ${TSTART} < traj.lammpstrj | ./raw2blockaverage.py 2 \
#             > traj_center-of-mass_trace.raw
#./raw2subtractlines.py -norm < traj_center-of-mass_trace.raw \
#             > traj_polymer_axis_directions.raw

# New method: Assume that the polymer axis is perpendicular to the direction
#             of the "bond" connecting two bases in opposite strands.  (I.E.
#             each "rung" of the "ladder" that represents the DNA helix.)
#             Hence we take the cross-product between successive "rungs" as
#             an approximation for the polymer axis direction at that location.
dump2data.py -raw -tstart ${TSTART} < traj.lammpstrj | ./merge_lines_periodic.py 0 1 -p 2 | awk '{if (NF==6) {x1=$1;y1=$2;z1=$3;x2=$4;y2=$5;z2=$6; print x2-x1" "y2-y1" "z2-z1} else print ""}' | ./merge_lines_periodic.py 0 1 -p 1 | awk '{if (NF==6) {x1=$1;y1=$2;z1=$3;x2=$4;y2=$5;z2=$6; x=y1*z2-y2*z1; y=z1*x2-z2*x1; z=x1*y2-x2*y1; r=sqrt(x*x+y*y+z*z); print x/r" "y/r" "z/r} else print ""}' > traj_polymer_axis_directions.raw
	    

# Download and the source code for the "ndautocrr_src" program 
# we will be using to calculate the correlation length of the polymer:

if [ ! -d "ndautocrr_src" ]; then
  git clone https://github.com/jewettaij/ndautocrr ndautocrr_src
fi

# Compile this code:
rm -f "ndautocrr"  # <- delete any old verions of this binary
cd ndautocrr_src/src/
  source setup_gcc.sh
  make clean
  make
  mv ndautocrr ../../
cd ../../

# Now run "ndautocrr" to calculate the auto-correlation function of the
# orientations of the polymer stored in "traj_polymer_axis_directions.raw".

./ndautocrr -L 16 -avezero \
            < traj_polymer_axis_directions.raw \
            > correlation_axis_orientation.dat


# The "ndautocrr" program above estimates the correlation length by estimating
# the area under the curve.  This is error-prone due to the large amount of
# noise at large separation lengths.
# For a more accurate estimate, examine the correlation function
# a short distance away, preferably at an integer number of periods (10.5bp)
# Try n=7  (Since 7*3 =21 ~= 2*10.5)
#     b=3*0.337 (where 3 is the number of base pairs per "monomers")
# We expect the correlation function to decay according to:
# C(n) = exp(-n/N)   (where N is the correlation length in units of monomers)
# solving for N = -i/log(C(i))

N=7
echo "#--------------------------------------------" >&2
echo "  more accurate estimate of persistence length (IN MONOMERS, not in nm):" >&2
awk "BEGIN{n=$N; b=1} {if ((NF==2) && (\$1==n)) {C_n=\$2; print b*(-n/log(C_n))}}" < correlation_axis_orientation.dat > stats_torsional_persistence_length_in_monomers_N=$N.dat
cat stats_torsional_persistence_length_in_monomers_N=$N.dat >&2
echo "#--------------------------------------------" >&2

