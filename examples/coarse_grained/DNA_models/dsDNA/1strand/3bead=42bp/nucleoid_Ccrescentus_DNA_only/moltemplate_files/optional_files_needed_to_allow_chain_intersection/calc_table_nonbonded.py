#!/usr/bin/env python

"""
Usage: 

  calc_table_nonbonded.py LABEL N rmin rmax eps sig Lcoeff rmaxLJ rshiftLJ Ldebye zi zj rshiftQ U0 \
       > table_DNA_backbone.dat

Example:

  ./calc_table_nonbonded.py BACKBONE 801 0.0 8.0 0.2 2.8 1 8.0 1.0 0.0 3 3 0.0 0.5961621 > table_DNA_backbone.dat
  ./calc_table_nonbonded.py BACKBONE 801 0.0 8.0 0.5961621 9.58 1 0  1.0 0.0 6 6 0.0 0.149040525 > table_DNA_backbone.dat

 Calculate a table of N pairwise energies and forces between 
 two charges qi=zi*qe, qj=zj*qe (qe=the magnitude of the charge of an electron)
 in a dielectric medium (assumed to be eps_r=80 for water), with Debye length L,
 and with an optional barrier height U0. (U0 is the maximum energy if two 
 particles lie on top of eachother.  It is infinity by default).

 This table will be read by LAMMPS using "pair_style table"
 http://lammps.sandia.gov/doc/pair_table.html

 Units are assumed to be in kcal/mole for energy
                       and nanometers for distance
    (and multiples of electron charge for charge)

 A constant is added to the potential energy so that it is 0 at rmax.
 Since LAMMPS extrapolates based on the last two entries,
 An extra entry in the table is added at the end with U=0 and F=0.
 This inures that both the force and the potential are zero at large r.
 This means the final table will have N+1 data points instead of N.


"""

comment_str = """
# The forces between two particles are the sum of Lennard-Jones
# interactions and Debye-screened coulombic interactions.
# Additionally, the potential is compressed so that the maximum 
# possible potential energy is finite (the \"U0\" parameter), and the 
# final potential is shifted so that the energy at rmax decays to
# zero at r = rmax
#
# Formula used:
#
#  By default, the energy between particles is defined as:
#      U(r) = U_LR(r) + U_yuk(r)
#
#  -- where --
#      U_LJ(r) is a generalized Lennard-Jones function:
#      U_LR(r) = ε * ((σ/(r-rshiftLJ))^12 - 2*λ*(σ/(r-rshiftLJ))^6)
#           Lennard-Jones parameters (ε, σ, λ, rshiftLJ),
#
#    (Special case: λ = 1,  rshiftLJ = 0
#     --> U_LJ(r) = 4ε((σ'/r)^12-(σ'/r)^6)    with σ' = σ*2^(-1/6) )
#
#  -- and --
#
#  U_yuk(r) = (Ke/r) * exp(-(r-rshiftQ)/Ldebye)     ("yukawa potential")
#           Ke = ke*qi*qj / eps_r
#           qi = the charge of the first particle
#           qj = the charge of the second particle
#        eps_r = is the relative dielectric permitivity != 80.0 for water
#       Ldebye = the Debye length   (The electrostatic decay length due to
#                counterions, such as electrolytes in the solution.  This is
#                set by the user.  In many biological systems it is ~=1nm)
#           ke = 8.9875517873681764e09 J*m*C^-2 = 1/4*pi*eps_0
#                (https://en.wikipedia.org/wiki/Coulomb's_law)
#                in units of (kCal/mole)*nm*e^-2
#      =8.9875517873681764e09*(6.0221409e+23/4184)*1e09*(1.6021766208e-19**2)/80
#                 Units: 1 J = (1/4184)*1kCal ("thermochemical calorie")
#                     1 kcal = 6.0221409e+23 kcal/mole
#                   and 1 nm = 1.0e-9 m
#                 (ke*e*e/eps_r = 0.415079644208)
#               (I should get something close to 0.71nm*kB*T, since
#               kB*T=0.001987207*300kCal/mole, and 0.71nm is the Bjerrum length)
#      rshiftQ = how close to the center are the charges in the object located?
#                (If on the surface, try rshiftQ = rshiftLJ + σ')
#
#  --- occlusion barrier: ---
#
#  If U0 the U0 parameter is not equal to float('inf'), then it calculates
#  this instead:
#            U(r) = (U0/(pi/2)) * arctan(Uorig(r)/U0)
#  where Uorig(r) = the original U(r) potential (shown above)
#              U0 = is the height of the energetic barrier for two particles
#                   to pass through each other
#              Note: This function approaches U0 as r->0
"""

import sys
from math import *

# derived constants:
qe    = 1.6021766208e-19   # the charge of an electron (in C)
eps_r = 80.0
ke    = 8.9875517873681764e09 # (= 1/4*pi*eps_0 =Coulomb's constant in J*m*C^-2)
Na    = 6.0221409e+23      # (Avogadro's number)
J_per_kcal = 4184          # J = (1/4184)*1kCal ("thermochemical calorie")
nm_per_m = 1.0e9
BIG_NUMBER = 1.0e100  # any large number will do

def U_yuk(r, sig, Ldebye, rshiftQ, zi, zj):
    Ke    = zi*zj*ke*(Na/J_per_kcal)*nm_per_m*(qe*qe)/eps_r
    # Ke = ke*qi*qj / eps_r   in units of (kcal/mole) * nm
    return (Ke/r) * exp(-(r-rshiftQ)/Ldebye)

def F_yuk(r, sig, Ldebye, rshiftQ, zi, zj):
    Ke    = zi*zj*ke*(Na/J_per_kcal)*nm_per_m*(qe*qe)/eps_r
    # Ke = ke*qi*qj / eps_r   in units of (kcal/mole) * nm
    return (Ke/r)*(1.0/r + (1.0/Ldebye)) * exp(-(r-rshiftQ)/Ldebye)

def U_LJ_unshifted(r, epsilon, sigma, Lcoeff):
    return epsilon * ((sigma/r)**12 - 2*Lcoeff*((sigma/r)**6))

def F_LJ_unshifted(r, epsilon, sigma, Lcoeff):
    return (epsilon/sigma) * (12*(sigma/r)**13 - 2*Lcoeff*6*((sigma/r)**7))

def U_LJ(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ):
    if (epsilon == 0):
        return 0.0
    elif (r < rshiftLJ):
        return BIG_NUMBER
    elif (r < rmaxLJ):
        return U_LJ_unshifted(r-rshiftLJ,      epsilon, sigma, Lcoeff) - \
               U_LJ_unshifted(rmaxLJ-rshiftLJ, epsilon, sigma, Lcoeff)
    else:
        return 0.0

def F_LJ(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ):
    if (epsilon == 0):
        return 0.0
    elif (r < rshiftLJ):
        return 0.0
    elif (r < rmaxLJ):
        return F_LJ_unshifted(r-rshiftLJ, epsilon, sigma, Lcoeff)
    else:
        return 0.0

def U_tot(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ):
    return U_LJ(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ) + \
           U_yuk(r, sigma, Ldebye, rshiftQ, zi, zj)

def F_tot(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ):
    return F_LJ(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ) + \
           F_yuk(r, sigma, Ldebye, rshiftQ, zi, zj)


def U(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ, U0, rmax):
    if (r > 0.0):
        uy = U_tot(r,   epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ)\
              - \
             U_tot(rmax,epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ)
        if U0 == 0.0 or uy < 0.0:
            return uy
        else:
            return (U0/(pi/2)) * atan(uy/(U0/(pi/2)))
    else:
        return U0


def F(r, epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ, U0, rmax):
    if (r > 0.0):
        uy = U_tot(r,   epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ) \
              - \
             U_tot(rmax,epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ)
        fy = F_tot(r,   epsilon, sigma, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ)
        if U0 == 0.0 or uy < 0.0:
            return fy
        else:
            u = uy / (U0/(pi/2))
            du_dr = -fy / (U0/(pi/2))
            datanu_du = 1.0/(1+u*u)    # derivative of arctan(u)
            datanu_dr = datanu_du * du_dr
            return -(U0/(pi/2)) * datanu_dr
    else:
        return 0.0


def PrintComments(out_file):
    out_file.write(comment_str)
    out_file.write("#\n"
                   "#\n"
                   "# columns:\n"
                   "#      i       r_i       U(r_i)     -dU/dr|r_i\n\n")


def PrintTable(out_file,
               eps, sig, Lcoeff, rmaxLJ,
               rshiftLJ, Ldebye, zi, zj, rshiftQ,
               U0,
               rmax, rmin, Ntable, label_str):

    if eps == float('inf'):
        eps = BIG_NUMBER
    if U0 == float('inf'):
        U0 = BIG_NUMBER

    #LAMMPS does not allow us to set rmin to 0
    if rmin == 0.0:
        rmin = rmax / 1.0e5

    if zi*zj == 0.0:
        rmax = rmaxLJ  #(use the smaller rmaxLJ cutoff)

    rminsq = rmin**2
    rmaxsq = rmax**2

    N = Ntable
    out_file.write(label_str+"\n")
    #rcut = rmin + N*(rmax-rmin)/(N-1)
    rcut = rmax
    out_file.write("N "+str(N)+" RSQ "+str(rmin)+" "+str(rcut)+#"\n")
                   " FPRIME 0.0 0.0\n\n")

    for i in range(0, N):
        rsq = rminsq + i*(rmaxsq-rminsq)/(N-1)
        r = sqrt(rsq)
        U_r = U(r, eps, sig, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ, U0, rmax)
        F_r = F(r, eps, sig, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ, U0, rmax)
        # Add an extra entry to the table at the end:
        if i == N:
            U_r = 0.0
            F_r = 0.0
        out_file.write(str(i+1)+' '+str(r)+' '+str(U_r)+' '+str(F_r)+'\n')
    out_file.write("\n")



if __name__ == "__main__":
    # argument list:
    #
    #  calc_table_nonbonded.py LABEL N rmin rmax eps sig Lcoeff rmaxLJ rshiftLJ Ldebye zi zj rshiftQ U0
    #

    label_str = sys.argv[1]
    N      = int(sys.argv[2])
    rmin   = float(sys.argv[3]) # typical value=0.0 (if U0 finite), 0.01 otherwise
    rmax   = float(sys.argv[4]) # typical value ~= 4.0 (assuming nm units)

    eps    = float(sys.argv[5])
    sig    = float(sys.argv[6])
    Lcoeff = 1.0
    if len(sys.argv) > 7:
        Lcoeff = float(sys.argv[7])
    rmaxLJ = rmax
    if len(sys.argv) > 8:
        rmaxLJ = float(sys.argv[8])

    rshiftLJ = 0.0
    if len(sys.argv) > 9:
        rshiftLJ = float(sys.argv[9])

    Ldebye = 1.0
    if len(sys.argv) > 10:
        Ldebye = float(sys.argv[10]) # Debye length (typically near 1.0nm)
    zi     = 0.0
    if len(sys.argv) > 11:
        zi=float(sys.argv[11])  # charge of particle i (in units of electron charge)
    zj     = 0.0
    if len(sys.argv) > 12:
        zj=float(sys.argv[12])  # charge of particle j (in units of electron charge)
    rshiftQ = 0.0
    if len(sys.argv) > 13:
        rshiftQ = float(sys.argv[13])
    U0     = 0.0
    if len(sys.argv) > 14:
        U0 = float(sys.argv[14]) # barrier height (for a single pair of particles)

    PrintTable(sys.stdout, eps, sig, Lcoeff, rmaxLJ, rshiftLJ, Ldebye, zi, zj, rshiftQ, U0, rmax, rmin, N, label_str)
