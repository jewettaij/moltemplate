# In the "run.in" file, I use these settings:
#
#                     kB*T kB*T tdamp  seed
fix fxlan all langevin 1.0 1.0  100.0 48279

# What "tdamp" parameter should we use?
#
# Recal "tdamp" is the time necessary for a particle to forget its former
# velocity due to random thermal/viscous forces (interactings with a heat bath).
# I sometimes complain about using large "tdamp" parameters.
# However when the polymer moves diffusively (in other words, when the 
# movement that we care about occur on timescales much larger than tdamp),
# THEN the system evolves at a rate which is inversely proportional to "tdamp".
# Consequently, for my initial tests I usually use a pretty large tdamp value.
# So we set the "tdamp" parameter to something similar to the save interval
# (used in the "dump" command below), multiplied by the "timestep" size,
# since there is no point in setting it to larger than this interval.
# The point here is not to reproduce physically accurate dynamics
# at short time intervals (which would require a very small tdamp
# parameter, like "1.0").  The point is to search conformational space
# efficiently at long, long timescapes.
# Later, perhaps we can reduce tdamp to 10.0 or 1.0, to see if it makes
# a difference (I can't imagine it does, but we can check this).
#
#
#
#
#
#
#
################################################################
# NOTES to myself explaining the parameters I (might be) using
################################################################
##################### DETAILS/SCRATCHWORK: ####################
#
# Choosing the parameters minimalist simulation of identical particles.
# which interact repulsively according to U(r)  ~  ε*(σ/r)^n 
#                                         ("ε"<-->"epsilon", "σ"<-->"sigma")
# Recall for diffusive motion:
#
#  <R^2> = D * t 
#      D = kB*T / gamma  (gamma = damping coeff w units force/velocity)
#        <==>
#  <R^2> = t * tdamp*(kB*T/m)        (tdamp = (gamma / m)^-1 = damping time)
#        <==>
#      D = tdamp*(kB*T/m)
#
#####
#
# Check: This makes sense in terms of a random walk.
#  <R^2> = N_s * b^2  (where N_s = number of steps, and b = the step size)
#        = (t / tdamp) * (v * tdamp)^2     ("v" = velocity, where m*v^2 ~ kB*T)
#        = t*  v^2  * tdamp
#        = t * (kB*T/m)*tdamp   (exploiting m*v^2 ~ kB*T)
#          (all this is ignoring dimensionless constants of order unity)
#
#####  Now introduce a simulation timestep "dt".  
#   Let t = N * dt
#
#  <R^2> = D * N * dt
#
#####  What is "dt"?
#
#   dt ~ w^(-1) = (k/m)^(-1/2)   ("w" = harmonic oscillation frequency
#                                 "k" = second derivative of U(r) around r)
#  Let U(r) = ε * (σ/r)^a
#
#        d^2          ε    / σ \^(a+2)
#    k = ---  U(r) = --- * |---|
#        dr^2        σ^2   \ r /
#
#  Where "r" is the distance of closest approach (between two particles).
#  At thermal equilibrium, this should satisfy
#
#  kB*T  ~  U(r)  =  ε * (σ/r)^a
#
#           / r \   /  ε \^(1/a)
#       --> |---| = |----|
#           \ σ /   \kB*T/
#
#                         [      m  / ε  \^(a+2)/a ] ^ (1/2)
# dt  ~  (m/k)^(1/2)  ~   [ σ^2 --- |----|         ]
#                         [      ε  \kB*T/         ]
#
##### Plug back into <R^2> = D N dt
# What we want to maximize is the distance traveled by diffusion
# (relative to the particle diameter, "σ")
#
#             _______D______  *  N  *  ________________dt_________________
#            '              `         '                                   `
#
# <R^2>    1   kB*T                     [ σ^2  m  / ε  \^(a+2)/a ] ^ (1/2)
# ----- = --- ------ * tdamp  *  N  *   [ ------- |----|         ]
#  σ^2    σ^2   m                       [    ε    \kB*T/         ]
#
# Constraints:  WLOG, set σ = 1
#                         ε = kB*T  (otherwise, "σ" is meaningless)
#
# Express <R^2>/σ^2 in terms of innertial time scale t_m = sqrt(m*σ^2 / kB*T)
#   Interpretation:
# t_m = time for a particle of kinetic-energy energy kB*T and mass m to travel s
#       (ballistically)
#
# <R^2>        1                                                  tdamp
# -----  = --------  *  tdamp  *  N  *  [ t_m^2 ] ^ (1/2)    =   --------  *  N
#  σ^2       t_m^2                                                 t_m
#
#
# If σ=1 and kB*T=1,
#
#        = (tdamp/sqrt(m)) * N   (subject to constraint that dt ~ sqrt(σ^2 m/ε))
#
# 
# RANGE CONSTRAINTS:
#
# 1)  dt must be << t_m  (the innertial time)
#     This happens autmatically because both t_m and dt scale  ~  sqrt(m)
#
# 2)  tdamp can not exceed "t"
#     (otherwise the simulation is not behaving diffusively on timescale t)
#
# So we set tdamp ~= t = N * dt
# This results in:
# <R^2>    1             1                 1  kB*T               1  kB*T
# ----- = ---  D * t  = --- D * N * dt  = --- ---- tdamp N*dt = --- ---- (N*dt)^2
#  σ^2    σ^2           σ^2               σ^2  m                σ^2  m 
#
#                                                                 1
#                                                             = ------  N^2 * dt^2
#                                                               t_m^2
#
#  However both t_m and dt are proportional to sqrt(m).
#
#  ---> Hence, the result (<R^2>) does not depend on the choice of m.
#       (assuming dt is chosen optimally).  This makes sense.
#       Under optimal conditions <R^2> only depends on the number
#       of timesteps, N, not m.
#  So, for convenience, choose m so that timescales are convenient.
#
###################################
# Recap:
#       σ  = in nm (varies, but near unity)
#       ε  = in kCal/mole (neary unity but varies depending on interaction)
#      kB  = 1  #(exactly, equivalent to "units lj", which we are using)
#       T  = 0.5961621 (0.001987207*300 in kCal/mole)
#      dt  = 1
#       m  = depends on dt and particle type (see below)
#       m  = (d^2 U(r)/dr^2) * (dt*Nperiod/2*pi)^2
#    tdamp ~= t  =  N * dt   # N~=simulation length
#        (solve for t_D.  Check that it should be << t)
#
#     In the past I used to set m=1 for convenience.
#     However, when multiple different kinds of particles and interactions are
#     present, I think it makes more sense to make dt a constant, and choose
#     the particle mass accordingly to maintain the same oscillation frequency.
#     (We are not trying to get microscopic dynamics correct, only the
#      thermodynamics, so we set all these frequencies to 1/Nperiod*dt)
#     Recall that:
#     dt = sqrt( m / (d^2 U(r)/dr^2) ) * (2*pi / Nperiod)
#      (Nperiod >> 1, dimensionless.  Typically, I set Nperiod = 10 or 20)
#      Check the force-field for the maximum value of (d^2U(r)/dr^2) 
#      which are accessible to the particle at the given temperature.
#      (For the specific case of U(r) ~= 1/r^n, then 
#       dt~=sqrt(σ*σ*m/kB*T) / Nperiod, but our case is more general)
#
#
#
# ------------ t_D (the diffusion timescale) ----------
# Useful to check what is t_D.
#        t_D = timescale over which diffusion displaces particles by 
#              a distance (squared) of σ^2. It is given by σ^2= D*t_D
#              (reminder: D = tdamp*(kB*T/m))
#  <-->  t_D = σ*σ / D
#
#        (The goal is to maximize <R^2>, not keep t_D = 1)
#
#        t_D should be << N*dt.   Make sure that this is so.  If not increase N
#
#        recall that 
#          D = (kB*T / m) * tdamp
