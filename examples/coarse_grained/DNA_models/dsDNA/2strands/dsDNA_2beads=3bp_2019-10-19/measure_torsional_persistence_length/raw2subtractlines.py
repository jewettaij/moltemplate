#!/usr/bin/env python


err_msg = """
Typical Usage:

   raw2subtractlines.py [-norm] < coordinate_file

   Coordinates read from one line of the file are subtracted from coordinates
   from the next line of the file (if it contains coordinates) and printed to
   the standard output.  Blank lines in the input file are copied to the
   standard out.  Each block of N lines of text containing M columns in the
   input file produces a block of N-1 lines of text (containing M columns)
   in the output file.

   The optional "-norm" argument allows you to normalize the resulting vectors
   after they have been subtracted.

Examples:

   raw2subtractlines.py       < coord_bead_chain.raw > coords_bond_vector.raw
   raw2subtractlines.py -norm < coord_bead_chain.raw > coords_bond_direction.raw

"""


import sys
from math import *
#import numpy as np


def ProcessStructure(x_id, normalize=False):
    D = len(x_id[0])
    N = len(x_id)
    for i in range(0, N-1):
        for d in range(0, D):
            x_diff = [x_id[i+1][d] - x_id[i][d] for d in range(0,D)]
        if (normalize):
            x_diff_len = 0.0
            for d in range(0, D):
                x_diff_len += x_diff[d] * x_diff[d]
            x_diff_len = sqrt(x_diff_len)
            for d in range(0, D):
                x_diff[d] /= x_diff_len

        sys.stdout.write(str(x_diff[0]))
        for d in range(1, D):
            sys.stdout.write(' ' + str(x_diff[d]))
        sys.stdout.write('\n')


# Parse the argument list:
if (len(sys.argv) > 2):
    sys.stderr.write("Error:\n\nTypical Usage:\n\n"+err_msg+"\n")
    exit(1)

if ((len(sys.argv) == 2) and
    ((sys.argv[1] == '-h') or
     (sys.argv[1] == '-?') or
     (sys.argv[1] == '--help'))):
    sys.stderr.write("Error:\n\nTypical Usage:\n\n"+err_msg+"\n")
    exit(1)

normalize = False

if (len(sys.argv) == 2):
    if ((sys.argv[1] == '-n') or
        (sys.argv[1] == '-norm') or
        (sys.argv[1] == '-normalize')):
        normalize = True
    else:
        sys.stderr.write("Error: Unrecognized command line argument:\n"
                         "       \""+sys.argv[1]+"\"\n")
        exit(1)


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
            ProcessStructure(x_id, normalize)
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
    x_d = list(map(float, tokens))
    x_id.append(x_d)

if len(x_id) > 0:
    ProcessStructure(x_id, normalize)

