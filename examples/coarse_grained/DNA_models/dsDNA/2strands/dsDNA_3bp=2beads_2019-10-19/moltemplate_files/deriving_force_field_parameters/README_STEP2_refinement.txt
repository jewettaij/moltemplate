# The procedure described in README_STEP1 generates a crude initial guess
# for the force-field parameters of the simulation a biomolecule.  After
# this guess, it is necessary to refine the choice of force field parameters.
# We did this iteratively by running simulations, measuring the distances
# and angles, adjusting the parameters, running more simulations, etc...
#
# At each iteration, a long simulation was run.
# After each simulation the distances and angles between the corresponding
# beads in the coarse grained model were measured (and modeled by a normal
# distribution "θ0_sim", "σsim").  Then the force-field parameters were
# adjusted (using the procedure described below), and the cycle was repeated.
# The goal is to iterate until:
#  θ0_sim=θ0_target  and  σsim = σtarget   (+/- some tolerance)

############### ITERATION PROCEDURE ####################
# "θ0" represents the target average "θ" (angle) parameter
# If the target θ0 is two degrees larger than the measured θ0
# (for example) then we increase the theta0 parameter by that amount.
# Average angle and distance parameters corrected using:

θ0_new = θ0_old + (θ0_target - θ0_sim)

# The normal relationship between the Hooke's law constant
# Utheta = (k/2)*(θ-θ0)^2 is:
# k_param = kB*T/σ^2
# (Careful.  LAMMPS uses this convention instead: K*(θ-θ0)^2)
# If the fluctuations in theta that we measured in our latest
# simulation were 10% larger (x1.10) larger than the target, 
# then we increase the k_param by the square of that ratio (1.10^2)

k_new = k_old * (σsim/σtarget)^2

##############################################################################
#   The procedure above works well for distances and angles, however we did not
# follow for dihedral angles. While this procedure generates DNA with the 
# correct shape, it does not generate DNA with the correct mechanical properties
# (ie. stiffness).  Apparently, DNA is crystallized under conditions which
# make it more straight than would be in liquid water under physiological
# conditions.  Consequently the DNA models generated this way are too stiff.
#   To deal with this, we chose dihedral angle stiffness parameters
# (and angle stiffness parameters) which reproduced the measured
# persistence length and torsional persistence length of double-stranded DNA.
# (50nm and 111nm, respectively).  This was usually done by hand, iteratively.
# The resulting measured mechanical properties were:
#                                                |  this model  | experimenal |
# ----------------------------------------------------------------------------
# persistence length:                            |  48.7101 nm  |     50 nm   |
# torsional persistence length:                  | 109.9752 nm  |    110 nm   |
# helical twist angle between monomers (3bp):    | 102.7797 deg |  see below  |
# -> helical twist angle between base-pairs:     |  34.2599 deg |   34.3 deg  |
# distance along the axis between monomers (3bp):|  0.98293 nm  |  see below  |
# -> distance along the axis between base-pairs: |  0.32764 nm  |   0.332 nm  |
##############################################################################
