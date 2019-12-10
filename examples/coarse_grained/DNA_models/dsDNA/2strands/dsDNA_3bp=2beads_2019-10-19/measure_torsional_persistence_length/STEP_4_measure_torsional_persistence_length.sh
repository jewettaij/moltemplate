#!/usr/bin/env bash

#dump2data.py -raw -tstart 1 < traj.lammpstrj | ./raw2blockaverage.py 2 \
#             > traj_center-of-mass_trace.raw
#
# "-tstart 1" skips the first frame.  
# The simulation might need more than one dump interval of time before 
# equilibrating.  In that case, use a larger -tstart value.

N=7   # separation between monomers used to estimate torsional persist length.
NPARTICLES_PER_MONOMER=2
NBP_PER_MONOMER=3

TSTART=10000000  # How many verlet iterations of simulation do we throw away?

dump2data.py -raw -tstart $TSTART < traj.lammpstrj > traj.raw


./merge_lines_periodic.py -p $NPARTICLES_PER_MONOMER  0 1 < traj.raw | awk '{if (NF==6) {dx=$4-$1; dy=$5-$2; dz=0.0; r=sqrt(dx*dx+dy*dy+dz*dz); print dx/r" "dy/r" "dz/r} else {print ""}}' > traj_basepair_axis.raw

awk 'BEGIN{pi=3.141592653589793} {if (NF==3) {angle=atan2($2,$1)*180.0/pi; print angle} else {print ""}}' < traj_basepair_axis.raw > traj_basepair_angle.dat

ANGLE_AT_REST=102.857142  # = NBP_PER_MONOMER * 360/10.5
                          #   guess can be within +/-90 degrees of truth

awk "BEGIN{preva=0; prevNF=0} {if (NF==0) {preva=0; prevNF=0; print ""} else {a=\$1; if (a<preva) {a=a+360.0} if (prevNF>0) {da=a-preva; if (da-$ANGLE_AT_REST<-180.0) {da=da+360} if (da-$ANGLE_AT_REST>180.0) {da=da-360} print da}} preva=a; prevNF=NF}" < traj_basepair_angle.dat > traj_basepair_delta_angle.dat


# Now obtain a more accurate estimate of the average helical pitch of the chain
# per monomer.  (This helical pitch angle will be weakly effected by the finite 
# temperature of the chain.  It might also be effected by the fact that we have
# confined the chain to a thin tube. Hopefully the tube will not effect it much)
ANGLE_AT_REST=`awk '{if (NF>0) {sum=sum+$1; n++}} END{printf("%.11g\n",sum/n)}' < traj_basepair_delta_angle.dat`

ANGLE_PER_BP=`awk "BEGIN{print $ANGLE_AT_REST/$NBP_PER_MONOMER}"`

echo "helical rotation per base pair (monomer=${NBP_PER_MONOMER}bp) = $ANGLE_PER_BP" >&2



#awk "{if (NF==0) {print ""} else {a=\$1-$ANGLE_AT_REST; print a}}" < traj_basepair_delta_angle.dat > traj_basepair_delta_angle_minus_expected.dat


awk "BEGIN{print 0.0} {if (NF==0) {print ""; print 0.0; sum=0.0} else {a=\$1-$ANGLE_AT_REST; sum+=a; print sum}}" < traj_basepair_delta_angle.dat > traj_basepair_angle_minus_expected.dat


# Calculate the difference in angle between monomer i and monomer i+N
# Sum over all pairs with this separation length.
# (Play around with different N to see what works best.  No need to check all N)
# Return the standard deviation of this number.  (The average should be zero.)

CALC_STDDEV_ANGLE_SRC=$(cat <<EOF
import sys
N = float(sys.argv[1])
from math import *
if sys.version[0] <= '2':
    import Queue
    q = Queue.Queue()
else:
    import queue
    q = queue.Queue()
sumsq = 0.0
n = 0
for line in sys.stdin:
    tokens = line.strip().split()
    if len(tokens) == 1:
        angle = float(tokens[0])
        q.put(angle)
        if q.qsize() == N+1:
            angle_prev = q.get()
            delta_angle = angle - angle_prev
            sumsq += delta_angle**2
            n += 1
    else:
        q.queue.clear()
print(str(sqrt(sumsq / n))+' '+str(n))
#print(str(sumsq / n)+' '+str(n))
EOF
)

echo "$CALC_STDDEV_ANGLE_SRC" > calc_stddev_angle.py

python calc_stddev_angle.py $N < traj_basepair_angle_minus_expected.dat \
      | awk '{print $1}' > stats_fluctuation_in_angle_in_degrees_N=$N.dat


rm -f calc_stddev_angle.py



awk '{print $1*(3.141592653589793/180.0)}' < stats_fluctuation_in_angle_in_degrees_N=$N.dat > stats_fluctuation_in_angle_in_radians_N=$N.dat

SIGMA_DELTA_ANGLE=`cat stats_fluctuation_in_angle_in_radians_N=$N.dat`


# For reasons I don't understand, the temperature of the simulation
# is often not what we requested it to be.  The force-field parameters
# were assigned assuming the simulation temperature would be 300Kelvin.
# The (torsional) persistence lengths are inversely temperature dependent.
# If the actual temperature was different, we need to rescale the
# temperature by this temperature difference.

T_SIM=`awk -v TSTART=$TSTART 'function isnum(x){return(x==x+0)} {if ((NF==8) && isnum($1) && isnum($2) && isnum($3) && isnum($4) && isnum($5) && isnum($6) && isnum($7) && isnum($8) && ($1>TSTART)) {sum+=$2; n+=1}} END{print sum/n}' < log.lammps`

T_TARGET=`echo "0.001987207*300" | bc`

SCALE_BY_TEMP=`echo "$T_SIM $T_TARGET" | awk '{print($1/$2)}'`




# Persistence length C:
# I am using the definition of "C" provided by Marko PhysRevE 2007, equation 7
# (That equation describes energy-per-unit-length of a twisted chain, not
#  total energy of a chain of length L.  So I multiply Marko's formula by L
#  to get the total twist energy of a chain of length L twisted (uniformly)
#  by angle "delta_angle".  Don't worry about the assumption that the chain
#  is being twisted uniformly.  You can can split the energy into components
#  due to uniform twisting which only depend on the difference between the
#  angle at the two ends of the chain, and local fluctuations about that 
#  uniform twist angle.  For details, see equation 7 of:
#  Rangle et. al,  J. Phys. Chem. B, (2008), Vol. 112, No. 42, p. 13359
# 
# U(delta_angle) = (1/2)*(kB*T*C/L)*(delta_angle)^2    <--definition of "C"
#              L = the length of the polymer (constrained not to bend)
# Probability(delta_angle) ~ exp(-(C/(2*L))*(delta_angle)^2)
#                          ~ exp(-(1/2)*(delta_angle/sigma)^2)
#  Since "sigma" = standard_deviation_in_delta_angle,
#  Solving for C:
#    C = L/(sigma**2)

CALC_C_SRC=$(cat <<EOF
from math import *
import sys
sigma = float(sys.argv[1])
L = float(sys.argv[2])
C = L/(sigma**2)
C *= $SCALE_BY_TEMP
print(C)
EOF
)

echo "$CALC_C_SRC" > calc_C.py

python calc_C.py $SIGMA_DELTA_ANGLE $N > stats_torsional_persistence_length_in_monomers_N=$N.dat

echo "fluctuations_in_angle(N=$N) = $SIGMA_DELTA_ANGLE" >&2


# Chain length statistics (used to calculate axial length per monomer)

awk 'BEGIN{prevNF=0; xp=0; yp=0; zp=0} {if (NF==0) {if (prevNF!=0) {print sqrt((xb-xa)^2+(yb-ya)^2+(zb-za)^2)} prevNF=0} else {if (prevNF==0) {xa=$1; ya=$2; za=$3} xb=$1; yb=$2; zb=$3} prevNF=NF}' < traj.raw > stats_chain_lengths.dat

CHAIN_LENGTH_NM=`awk 'BEGIN{s=0;n=0} {if (NF!=0){s+=$1; n++}} END{print s/n}'<stats_chain_lengths.dat`

CHAIN_LENGTH_MONOMERS=`awk "BEGIN{first=1} {if ((NF==0) && (first)) {print (NR-1)/$NPARTICLES_PER_MONOMER; first=0}}" < traj.raw`

CHAIN_LENGTH_NM_PER_MONOMER=`awk "BEGIN{print $CHAIN_LENGTH_NM/$CHAIN_LENGTH_MONOMERS}"`

CHAIN_LENGTH_NM_PER_N_MONOMERS=`awk "BEGIN{print $N*$CHAIN_LENGTH_NM_PER_MONOMER}"`

CHAIN_LENGTH_NM_PER_BP=`awk "BEGIN{print $CHAIN_LENGTH_NM_PER_MONOMER/$NBP_PER_MONOMER}"`

echo "chain length per base pair (monomer=$NBP_PER_MONOMER) = $CHAIN_LENGTH_NM_PER_BP" >&2

python calc_C.py $SIGMA_DELTA_ANGLE $CHAIN_LENGTH_NM_PER_N_MONOMERS > stats_torsional_persistence_length_in_nm_N=$N.dat

rm -f calc_C.py


# As a sanity check, calculate the fluctuations in angle between successive 
# "monomers" (by setting L=b), and then compare this with the fluctuations 
# in successive "basepair" dihedral angles. (Again, each "basepair" in the model
# corresponds to several base-pairs in real double-stranded DNA; usually 3.)

#echo ' Angle fluctuations between successive \"monomers\":'
#awk '{C=$1; L=1.0108; kB=0.001987207; T=300; print sqrt(L/C)*180/3.141592653589793}' < stats_torsional_persistence_length_in_nm_N=$N.dat

echo ' The torsional persistence length of the polymer in nm is stored in this file:'
echo '    stats_torsional_persistence_length_in_nm_N=$N.dat'
cat stats_torsional_persistence_length_in_nm_N=$N.dat

