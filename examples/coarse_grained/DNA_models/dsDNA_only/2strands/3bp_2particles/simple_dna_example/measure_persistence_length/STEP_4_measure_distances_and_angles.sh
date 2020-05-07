#!/usr/bin/env bash


# We need to make sure "dlpdb" is installed.  (It has some useful commands
# for extracting distances and angles from multicolumn numeric text files.)
# If not, get it from pypi:

if ! type coords2distances.py > /dev/null; then
  DOWNLOADED_DLPDB="true"
  echo "installing dlpdb" >&2
  python -m venv venv
  source venv/bin/activate
  pip install dlpdb
fi



TSTART=10000000  # How many verlet iterations of simulation do we throw away?

dump2data.py -raw -tstart $TSTART < traj.lammpstrj > traj_t=${TSTART}-.raw

./merge_lines_periodic.py -p 2  0 2 < traj_t=${TSTART}-.raw | awk '{if (NF==6) print $0}' | coords2distances.py > stats_distances_Backbone_.dat
./merge_lines_periodic.py -p 2  1 3 < traj_t=${TSTART}-.raw | awk '{if (NF==6) print $0}'| coords2distances.py >> stats_distances_Backbone_.dat
./calc_stddev.awk stats_distances_Backbone_.dat > stats_distances_Backbone_ave_stddev.dat

./merge_lines_periodic.py -p 2  0 1 < traj_t=${TSTART}-.raw | awk '{if (NF==6) print $0}' | coords2distances.py > stats_distances_BasePair_.dat
./calc_stddev.awk stats_distances_BasePair_.dat > stats_distances_BasePair_ave_stddev.dat

# acute angles pointing in the forward direction (diction of increasing index)
./merge_lines_periodic.py -p 2  0 1 3 < traj_t=${TSTART}-.raw | awk '{if (NF==9) print $0}' | coords2angles.py > stats_angles_BasePairAcute_.dat
# ############
# acute angles pointing in the reverse direction (diction of decreasing index)
# COMMENTING OUT:
#./merge_lines_periodic.py -p 2  0 2 3 < traj_t=${TSTART}-.raw | awk '{if (NF==9) print $0}' | coords2angles.py >> stats_angles_BasePairAcute_.dat
# ############
# WHY:
# ############
# Ideally a model of double-stranded DNA (dsDNA) would be symmetric.  It would
# that look identical when viewed down the axis from either direction.  (dsDNA
# is symmetric because two DNA strands in dsDNA point in opposite directions.)
# In this model we are using a (slightly) asymmetric force-field.
# When the force-field not symmetric, it means that in principle a user might
# be able to detect that the DNA polymer looks slightly different viewing it
# along one direction compared to the reverse direction.  Why did we do this?
# The various angle and dihedral constraints overlap and are somewhat redundant.
# So in this particular model, we decided to omit some of them
# in order to reduce the computational cost of the force calculation.
# The force field we are using only directly applies a force to some of these
# acute angles (the forward acute angles) and not others (reverse acute angles).
# Since we use the results of these measurements to tune those force-field
# parameters, we prefer to only measure those acute angles we can control.
# So, I commented out the code that measures the reverse acute angles above.
# Now, calculate the standard deviation of these angle measurements.
./calc_stddev.awk stats_angles_BasePairAcute_.dat > stats_angles_BasePairAcute_ave_stddev.dat


./merge_lines_periodic.py -p 2  1 0 2 < traj_t=${TSTART}-.raw | awk '{if (NF==9) print $0}' | coords2angles.py > stats_angles_BasePairObtuse_.dat
./merge_lines_periodic.py -p 2  1 3 2 < traj_t=${TSTART}-.raw | awk '{if (NF==9) print $0}' | coords2angles.py >> stats_angles_BasePairObtuse_.dat
./calc_stddev.awk stats_angles_BasePairObtuse_.dat > stats_angles_BasePairObtuse_ave_stddev.dat


./merge_lines_periodic.py -p 2  1 0 2 3 < traj_t=${TSTART}-.raw | awk '{if (NF==12) print $0}' | coords2dihedrals.py 225.0 > stats_dihedrals_BasePair_.dat
./merge_lines_periodic.py -p 2  0 1 3 2 < traj_t=${TSTART}-.raw | awk '{if (NF==12) print $0}' | coords2dihedrals.py 225.0 >> stats_dihedrals_BasePair_.dat
./calc_stddev.awk stats_dihedrals_BasePair_.dat > stats_dihedrals_BasePair_ave_stddev.dat


./merge_lines_periodic.py -p 2  1 3 2 4 < traj_t=${TSTART}-.raw | awk '{if (NF==12) print $0}' | coords2dihedrals.py 225.0 > stats_dihedrals_MajorGroove_.dat
./calc_stddev.awk stats_dihedrals_MajorGroove_.dat > stats_dihedrals_MajorGroove_ave_stddev.dat

./merge_lines_periodic.py -p 2  2 0 1 3 < traj_t=${TSTART}-.raw | awk '{if (NF==12) print $0}' | coords2dihedrals.py 60 > stats_dihedrals_torsion_.dat
./calc_stddev.awk stats_dihedrals_torsion_.dat > stats_dihedrals_torsion_ave_stddev.dat



### CLEANUP (uninstall dlpdb)

if [ -n "$DOWNLOADED_DLPDB" ]; then
  echo "uninstalling dlpdb" >&2
  deactivate
  rm -rf venv
fi

