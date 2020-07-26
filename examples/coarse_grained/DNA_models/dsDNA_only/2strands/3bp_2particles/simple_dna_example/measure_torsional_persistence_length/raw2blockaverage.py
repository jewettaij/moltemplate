#!/usr/bin/env python


err_msg = """
Typical Usage:

   raw2blockaverage.py N [scale_inv] < coordinate_file

   Coordinates read from the file coordinate_file are averaged in blocks
   of size N, and printed to the standard output, followed by a blank line.
   Excluding blank lines, the number of lines in the output equals the number 
   of lines in the input divided by N.  If blank lines are present, then
   the coordinates read from the file are assumed to represent independent
   snapshots from a trajectory (animation).  In this case, the block-averaging
   is done repeatedly for each frame in the animation, and a new trajectory
   file is written (containing blank line delimters between frames).

   The optional "scale_inv" argument allows you to divide the 
   all of resulting averaged coordinates by the number scale_inv.
   (Typically, N and scale_inv, if present, are equal to each other.)

Example:

   raw2blockaverage.py 2 < coords.raw > coords_ave2.raw
   raw2blockaverage.py 3 3 < coords.raw > coords_ave3_normalized.raw

"""


import sys
from math import *
#import numpy as np


def ProcessStructure(x_id, n_ave, scale):
    D = len(x_id[0])
    n_orig = len(x_id)
    for i in range(0, n_orig/n_ave):
        xave_d = [0.0 for d in range(0, D)]
        for j in range(0, n_ave):
            for d in range(0, D):
                xave_d[d] += x_id[n_ave*i + j][d]
        for d in range(0, D):
            xave_d[d] *= scale/float(n_ave)
        sys.stdout.write(str(xave_d[0]))
        for d in range(1, D):
            sys.stdout.write(' '+str(xave_d[d]))
        sys.stdout.write('\n')


# Parse the argument list:
if len(sys.argv) <= 1:
    sys.stderr.write("Error:\n\nTypical Usage:\n\n"+err_msg+"\n")
    exit(1)

n_ave = int(sys.argv[1])

scale = 1.0
if len(sys.argv) > 2:
    scale = 1.0 / float(sys.argv[2])


# Now read the input file:

x_id = []
count_structs = 0
is_new_structure = True
interpret_blank_lines_as_new_structures = True

in_file = sys.stdin
for line_orig in in_file:

    ic = line_orig.find('#')
    if ic != -1:
        line = line_orig[:ic]
    else:
        line = line_orig.rstrip('\n')    
    
    tokens = line.strip().split()
    if len(tokens) == 0:
        if (interpret_blank_lines_as_new_structures and
            (len(x_id) > 0)):
            # blank (or comment) lines signal the next frame of animation
            ProcessStructure(x_id, n_ave, scale)
            sys.stdout.write('\n')

            x_id = []
            count_structs += 1
            #sys.stderr.write('done\n')
            is_new_structure = True
        continue  # skip blank lines or comments

    elif is_new_structure:
        is_new_structure = False

    # x_d contains the coordinates read from the 
    # most recent line in the current frame
    x_d = map(float, tokens)
    x_id.append(x_d)

if len(x_id) > 0:
    ProcessStructure(x_id, n_ave, scale)

