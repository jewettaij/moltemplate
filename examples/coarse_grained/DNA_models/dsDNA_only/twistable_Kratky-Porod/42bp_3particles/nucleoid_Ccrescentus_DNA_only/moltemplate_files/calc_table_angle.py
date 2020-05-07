#!/usr/bin/env python

# We keep track of the program name and version.
# (This is only used for generating error messages.)
g_program_name = __file__.split('/')[-1]
g_date_str = '2017-3-24'
g_version_str = '0.1.0'


import sys
from math import *

theta0      = float(sys.argv[1]) # position of minima
K           = float(sys.argv[2]) # second derivative at theta0
theta_range = float(sys.argv[3]) # energy->infinity as theta-theta0->theta_range
N           = int(sys.argv[4])   # number of entries in the table. typically 181
LABEL       = sys.argv[5]        # A name to give this table
theta_min = 0.0
if len(sys.argv) >= 6:
    theta_min   = float(sys.argv[6])
if len(sys.argv) >= 7:
    theta_max   = float(sys.argv[7])
theta_restrict_min = theta_min
F_restrict_min = 10.0
if len(sys.argv) > 9:
    theta_restrict_min = float(sys.argv[8])
    F_restrict_min = float(sys.argv[9])
theta_restrict_max = theta_max
F_restrict_max = -F_restrict_min
if len(sys.argv) > 11:
    theta_restrict_max = float(sys.argv[10])
    F_restrict_max = float(sys.argv[11])
OTHER_PARAM = ""                # other parameters following LABEL and N. See:
if len(sys.argv) >= 12:         # http://lammps.sandia.gov/doc/angle_table.html
    OTHER_PARAM = sys.argv[12]

# Calculate the energy (and energy derivative) of the bond-angle-potential
#
#   U(theta) = ((2/pi)t_range)^2 * K (-1+1/ cos((theta-theta0)/((2/pi)t_range)))
#
# http://lammps.sandia.gov/doc/angle_table.html
#
# Note: The second derivative at the minima is "K" (in radians)
#       In the absence of other forces, and in the limit of small K, this 
#       is also the persistence length (in units of backbone-bond-lengths).
#       But this is only true for simple linear (1-bead) chains with only
#       one particle per monomer, and no steric (size-exclusion) forces.


# The third column contains the energy as a function of theta (in degrees)
def U(theta, theta0, K, theta_range=180.0):
    #   U(theta)=((2/pi)t_range)^2*K*(-1+1/ cos((theta-theta0)/((2/pi)t_range)))
    t  = theta*pi/180.0
    t0 = theta0*pi/180.0
    t_range = theta_range*pi/180.0
    kappa = (2/pi)*t_range
    return (kappa**2) * K * ((1.0 / cos((t-t0)/kappa)) - 1.0)


# The fourth column contains the negative derivative of the energy as a function
# of theta. (Both theta and the derivative respect to theta are in degrees.)
def F(theta, theta0, K, theta_range=180.0):
    #   U(theta) = ((2/pi)t_range)^2 * K / cos((theta-theta0)/((2/pi)t_range))
    t  = theta*pi/180.0
    t0 = theta0*pi/180.0
    t_range = theta_range*pi/180.0
    kappa = (2/pi)*t_range
    return -kappa * K * sin((t-t0)/kappa) / ((cos((t-t0)/kappa))**2)


sys.stdout.write('# This file was created using:\n'
                 '# '+g_program_name+' '+' '.join(sys.argv[1:])+'\n'
                 '#\n')

sys.stdout.write("""
#   U(theta)=((2/pi)t_range)^2 * K * (-1+1/ cos((theta-theta0)/((2/pi)t_range)))
#
""")


sys.stdout.write('#    theta0      = ' + str((pi/180.0)*theta0) + ' in radians  (' + 
                 str(theta0) + ' in degrees)\n' +
                 '#    t_range = ' + str((pi/180.0)*theta_range) + ' in radians  (' + 
                 str(theta_range) + ' in degrees)\n'+
                 '#    K         = ' + str(K) + ' (in energy/radian*2)\n'
                 '#\n'
                 '# Table format:\n'
                 '#\n'
                 '# i  theta_i  U(theta_i)  -dU/dtheta|theta_i   (in degrees)\n'
                 '\n' +
                 LABEL + '\n'
                 'N '+str(N)+' '+OTHER_PARAM+'\n'
                 '\n')

for i in range(0,N):
    theta = theta_min + i*(theta_max - theta_min)/(N-1)
    if theta < theta_restrict_min:
        f = F_restrict_min*180/pi
        u = U(theta_restrict_min, theta0, K, theta_range) - \
            (theta-theta_restrict_min)*f*pi/180
    elif theta > theta_restrict_max:
        f = F_restrict_max*180/pi
        u = U(theta_restrict_max, theta0, K, theta_range) - \
            (theta-theta_restrict_max)*f*pi/180
    else:
        u = U(theta, theta0, K, theta_range)
        f = F(theta, theta0, K, theta_range)
    sys.stdout.write(str(i+1)+' '+str(theta)+' '+str(u)+' '+str(f*pi/180)+'\n')


