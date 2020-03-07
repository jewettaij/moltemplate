#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

g_program_name = __file__.split('/')[-1]  # = 'interpolate_curve.py'
g_version_str = '0.3.0'
g_date_str = '2019-12-12'

g_usage_str = """
Usage:

   """ + g_program_name + """ Ndesired [scale] [alpha] < coords_orig.raw > coords.raw

Example:

   """ + g_program_name + """ 30117 3.0 0.5 < coords_orig.raw > coords.raw

"""


import sys
from math import *
import numpy as np
import bisect



## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver stolen from:
# https://gist.github.com/cbellei/8ab3ab8551b8dfc8b081c518ccd9ada9
# (also taken from https://gist.github.com/ofan666/1875903)

import numpy as np

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    Solves a (non-cyclic) system of equations of the form:

      a_i*x_{i-1} + b_i*x_i + c_i*x_{i+1} = d_i
      where a_1=0, and c_n=0

    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc



def CalcNaturalCubicSplineCoeffs(r, spline_exponent_alpha=0.5):
    N = len(r)
    assert(N >= 4)
    D = len(r[0])
    tcontrol = np.zeros(N)

    # e_d2rdt2[d][i] is the second derivative of the spline at the ith 
    #               control point (and in the dth direction)
    #
    # Once we have figured out e_d2rdt2[d][i], we can calculate the spline
    # anywhere using:
    #   SplineEval(t, h_i, 
    #              e_{i+1} / (6*h_i), 
    #              e_i / (6*h_{i+1}),
    #              r_{i+1}/h_i - e_{i+1}*h_i,
    #              r_i/h_i - e_i*h_i)

    e_d2rdt2 = np.zeros((D,N))

    # We want to solve this system of equations for e_1, e_2, ..., e_{n-2}:
    #    (note: e_i is shorthand for e_d2rdt2[d][i])
    #
    # h_{i-1}*e_{i-1} + u_i*e_i + h_{i+1} * e_{i+1}  =  v_i
    #   where h_i, u_i and v_i are shorthand for:
    #
    # h_i = change in "time" parameter for the ith interval, which is usually:
    #     = |r_{i+1} - r_i |^alpha    (alpha varies depending on settings)
    #   and
    # u_i = 2 * (h_{i-1} + h_i )
    # v_i = 6 * (b_i     - b_{i-1})
    # b_i =  (r_i - r_{i-1}) / h_i
    #
    # ...subject to the constraints that the curve is straight at the endpoints:
    # e_0 = 0       <-- first control point  (indexing begins at 0)
    # e_{n-1} = 0   <-- this is the last control point  (indexing begins at 0)

    h_dt = np.zeros(N-1)   # h_dt[i] is the i'th time interval 
                         # in the parameterization
    b_drdt = np.zeros((N-1,D))  # b_drdt is a discrete version of the derivative

    t_total=0.0
    for i in range(0, N-1):
        tcontrol[i] = t_total
        sqr_distance_i_ip1 = 0.0
        for d in range(0, D):
            sqr_distance_i_ip1 += (r[i+1][d] - r[i][d])**2
        h_dt[i] = sqr_distance_i_ip1**(0.5*spline_exponent_alpha)
        for d in range(0, D):
            b_drdt[i][d] = (r[i+1][d] - r[i][d]) / h_dt[i]
        t_total += h_dt[i]
    tcontrol[N-1] = t_total

    a_coeff = np.zeros(N)
    b_coeff = np.zeros(N)
    c_coeff = np.zeros(N)
    d_coeff = np.zeros((D,N))

    for i in range(1, N-1):
        # h_dt[i] is the difference in "time" in the parametric curve 
        # between pairs of control points.  If spline_exponenent_alpha is 0
        # then the time interval between control points is uniform.
        # ("spline_exponent_alpha" is 0.5 for centripital Catmull-Rom splines.)
        a_coeff[i]    =      h_dt[i-1]
        b_coeff[i]    = 2.0*(h_dt[i-1] + h_dt[i])
        c_coeff[i]    =      h_dt[i]
        for d in range(0, D):
            d_coeff[d][i] = 6.0*(b_drdt[i][d] - b_drdt[i-1][d])

    a_coeff[0] = 0.0
    b_coeff[0] = 1.0
    c_coeff[0] = 0.0
    for d in range(0, D):
        d_coeff[d][0] = 0.0
    a_coeff[N-1] = 0.0
    b_coeff[N-1] = 1.0
    c_coeff[N-1] = 0.0
    for d in range(0, D):
        d_coeff[d][N-1] = 0.0

    for d in range(0, D):
        e_d2rdt2[d] = TDMAsolver(a_coeff, b_coeff, c_coeff, d_coeff[d])

    # alternately, if that fails, try the matrix inverter that comes with numpy:
    #M = np.zeros((N,N))
    #for i in range(0,N):
    #    if i-1>=0:
    #        M[i][i-1] = a_coeff[i]
    #    M[i][i]   = b_coeff[i]
    #    if i+1 < N:
    #        M[i][i+1] = c_coeff[i]
    #for d in range(0, D):
    #    e_d2rdt2[d] = np.linalg.solve(M, d_coeff[d])

    e_d2rdt2 = e_d2rdt2.transpose()
    c3a = np.zeros((N, D))
    c3b = np.zeros((N, D))
    c1a = np.zeros((N, D))
    c1b = np.zeros((N, D))
    # faster to precompute these coefficients in advance:
    for i in range(0, N-1):
        # c3a = e_{i+1} / (6*h_i)
        # c3b =     e_i / (6*h_{i+1})
        # c1a = r_{i+1}/h_i - e_{i+1}*h_i
        # c1b =  r_i / h_i - e_i*h_i)
        c3a[i] = e_d2rdt2[i+1] / (6*h_dt[i])
        c3b[i] = e_d2rdt2[i]   / (6*h_dt[i])
        c1a[i] = r[i+1]/h_dt[i]  -  e_d2rdt2[i+1]*h_dt[i]/6.0
        c1b[i] = r[i]/h_dt[i]    -  e_d2rdt2[i]*h_dt[i]/6.0

    # Return these spline coefficients to the caller.
    # We can use these to quickly evaluate the spline repetatively later on

    return c3a, c3b, c1a, c1b, tcontrol





def SplineEval(t, t_interval, c3a, c3b, c1a, c1b):
    ta = t
    tb = t_interval - t
    return (c3a*ta*ta + c1a)*ta + (c3b*tb*tb + c1b)*tb



def SplineEvalD1(t, t_interval, c3a, c3b, c1a, c1b):
    ta = t
    tb = t_interval - t
    return c3a*ta*ta*3.0 + c1a - (c3b*tb*tb*3.0 + c1b)



def SplineEvalD2(t, t_interval, c3a, c3b, c1a, c1b):
    ta = t
    tb = t_interval - t
    return c3a*ta*6.0 + c3b*tb*6.0



def SplineInterpEval(t, c3a, c3b, c1a, c1b, tcontrol):
    i = bisect.bisect(tcontrol, t) - 1
    if i < 0: i = 0
    if i >= len(tcontrol)-1: i = len(tcontrol)-2
    return SplineEval(t-tcontrol[i], tcontrol[i+1]-tcontrol[i],
                      c3a[i], c3b[i], c1a[i], c1b[i])


def SplineInterpEvalD1(t, c3a, c3b, c1a, c1b, tcontrol):
    i = bisect.bisect(tcontrol, t) - 1
    if i < 0: i = 0
    if i >= len(tcontrol)-1: i = len(tcontrol)-2
    return SplineEvalD1(t-tcontrol[i], tcontrol[i+1]-tcontrol[i],
                        c3a[i], c3b[i], c1a[i], c1b[i])



def SplineInterpEvalD2(t, c3a, c3b, c1a, c1b, tcontrol):
    i = bisect.bisect(tcontrol, t) - 1
    if i < 0: i = 0
    if i >= len(tcontrol)-1: i = len(tcontrol)-2
    return SplineEvalD2(t-tcontrol[i], tcontrol[i+1]-tcontrol[i],
                        c3a[i], c3b[i], c1a[i], c1b[i])



def SplineCurvature2D(t, t_interval, c3a, c3b, c1a, c1b):
    # first derivatives
    x1,y1 = SplineEvalD1(t, t_interval, c3a, c3b, c1a, c1b)
    # second derivatives
    x2,y2 = SplineEvalD2(t, t_interval, c3a, c3b, c1a, c1b)
    denom = pow((x1*x1 + y1*y1), 1.5)
    curvature = (x1*y2 - x2*y1) / denom
    return curvature



def SplineInterpCurvature2D(t, c3a, c3b, c1a, c1b, tcontrol):
    i = bisect.bisect(tcontrol, t) - 1
    if i < 0: i = 0
    if i >= len(tcontrol)-1: i = len(tcontrol)-2
    return SplineCurvature2D(t-tcontrol[i], tcontrol[i+1]-tcontrol[i],
                             c3a[i], c3b[i], c1a[i], c1b[i])



def ResampleCurve(x_orig, num_points, alpha=0.5):
    """
       Given a list (or numpy array) of points in n-dimensional space that lie
    along some curve, this function returns a new list of "num_points" points
    which are uniformly distributed along the original curve.  This is done
    using cubic spline interpolation (which is tuned by the "alpha" parameter).

    https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull%E2%80%93Rom_spline

       Here "uniformly" means distributed at even intervals in the spline-
    parameter-space, not in physical distance. (If the original points were
    not uniformly distributed along the curve, then new points won't be either.)

       This function has not been optimized for speed.
    """

    assert(len(x_orig) >= 4)
    num_dimensions = len(x_orig[0])
    x_new = np.zeros((num_points, num_dimensions))

    c3a, c3b, c1a, c1b, tcontrol = \
        CalcNaturalCubicSplineCoeffs(x_orig, alpha)
    tmin = 0.0
    tmax = tcontrol[-1]
    for i in range(0, num_points):
        t = tmin + i*((tmax - tmin) / (num_points-1))
        x_new[i] = SplineInterpEval(t, c3a, c3b, c1a, c1b, tcontrol)

    return x_new


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return self.err_msg

    def __repr__(self):
        return str(self)


def main():

    try:

        #######  Main Code Below: #######
        sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

        spline_exponent_alpha = 0.5
        use_linear_interpolation = False

        # Parse the argument list:
        if len(sys.argv) <= 1:
            raise InputError("Missing arguments\n"+g_usage_str+"\n")
                             

        n_new = int(sys.argv[1])

        if len(sys.argv) > 2:
            scale = float(sys.argv[2])
        else:
            scale = 1.0

        if len(sys.argv) > 3:
            spline_exponent_alpha = float(sys.argv[3])            

        coords_orig = []

        lines = sys.stdin.readlines()

        for line in lines:
            tokens = line.split()
            if (len(tokens) > 0):
                coords_orig.append(list(map(float, tokens)))
                g_dim = len(tokens)

        n_orig = len(coords_orig)

        x_orig = np.array(coords_orig)

        if n_orig < 2:
            raise InputError("Input file contains less than two lines of coordinates.")

        if n_orig < 4:
            use_linear_interpolation = True

        if n_new < 2:
            raise InputError("Output file will contain less than two lines of coordinates.")

        #coords_new = [[0.0 for d in range(0, g_dim)] for i in range(0, n_new)]
        x_new = np.zeros((n_new, g_dim))

        if use_linear_interpolation:
            for i_new in range(0, n_new):
                I_orig = (i_new) * (float(n_orig-1) / float(n_new-1))
                i_orig = int(floor(I_orig))
                i_remainder = I_orig - i_orig

                if (i_new < n_new-1):
                    for d in range(0, g_dim):
                        x_new[i_new][d] = scale*(x_orig[i_orig][d]
                                                 +
                                                 i_remainder*(x_orig[i_orig+1][d]-
                                                              x_orig[i_orig][d]))
                else:
                    for d in range(0, g_dim):
                        x_new[i_new][d] = scale*x_orig[n_orig-1][d]

        else:
            x_new = ResampleCurve(x_orig, n_new, spline_exponent_alpha)

        # print the coordates
        for i in range(0, n_new):
            for d in range(0, g_dim-1):
                sys.stdout.write(str(x_new[i][d]*scale) + ' ')
            sys.stdout.write(str(x_new[i][g_dim-1]*scale) + "\n")

    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + 'Error:\n\n' + str(err) + '\n')
        sys.exit(1)


if __name__ == '__main__':
    main()
