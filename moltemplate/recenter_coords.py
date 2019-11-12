#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2015


g_program_name = __file__.split('/')[-1]  # "recenter_coords.py"
g_version_str  = '0.1.0'
g_date_str     = '2019-8-22'



g_usage_str = \
"""

   This program reads a (typically 3-column) ASCII coordinate file 
   (.raw format) and creates a new coordinate file where the average position
   (geometric center) is located at a position specified by the user.
   This program works in an arbitrary number of dimensions.
   The number of arguments must equal the number of columns in the input file.
   (This program does not understand blank lines.
    It does not work with blank-line-delimited trajectories.)

Usage:

   """ + g_program_name + """ Xcen [Ycen Zcen ...] < coords_old.raw > coords_new.raw

Example Usage: 

   """ + g_program_name + """ 0 0 0 < coords_old.raw > coords_new.raw

"""


from math import sqrt
import sys



class InputError(Exception):
    """ A generic exception object containing a string for error reporting.

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)

# x_id        stores the contents of the file (numeric multidimensional table)
#              (It is a 2-D table.)
# x_id[i]     stores a list of coordinates for the ith particle in the polymer
# x_id[i][d]  stores the d'th coordinate for the ith particle
#              (The number of coordinates need not be 3.)

x_id=[]

def main():
    try:
        #######  Main Code Below: #######
        sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

        if len(sys.argv) <= 1:
            raise(InputError('Error: Expected at least 1 argument\n'+g_usage_str))

        # Convert the arguments to numbers
        x_cm_new_d = list(map(float, sys.argv[1:]))
        g_dim = len(x_cm_new_d)

        for line_orig in sys.stdin:
            ic = line_orig.find('#')
            if ic != -1:
                line = line_orig[:ic]
            else:
                line = line_orig.rstrip('\n')    

            tokens = line.strip().split()
            if len(tokens) == 0:
                continue  # skip blank lines or comments

            D = len(tokens) # = number of coordinates on every line of the file

            x_d = list(map(float, tokens))
            x_id.append(x_d)

        N = len(x_id)
        sys.stderr.write('(done reading file)\n')

        # calculate the geometric center (average position)
        xcm_d = [0.0 for d in range(0,g_dim)]
        for i in range(0, N):
            for d in range(0, g_dim):
                xcm_d[d] += x_id[i][d]

        for d in range(0, g_dim):
            xcm_d[d] /= len(x_id)

        for i in range(0, N):
            sys.stdout.write(str(x_id[i][0] - xcm_d[0]))
            for d in range(1, g_dim):
                sys.stdout.write(" "+str(x_id[i][d] + x_cm_new_d[d] - xcm_d[d]))
            sys.stdout.write("\n")

    except (ValueError, InputError) as err:
        sys.stderr.write('\n\n'+str(err)+'\n')
        sys.exit(-1)



if __name__ == '__main__':
    main()
