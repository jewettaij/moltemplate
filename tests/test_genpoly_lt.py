#!/usr/bin/env python3

import math
import numpy as np
import moltemplate

N = 4
# Generate a zig-zag curve containing N points
x_orig = np.array([[i, 0.5*(i%2), 0.0] for i in range(0,N)])

# It's a really good idea to generate a smoother version of this curve:
x_new = moltemplate.interpolate_curve.ResampleCurve(x_orig, 21, 0.5)

# We want the spacing between monomers to be 0.332
x_new *= 0.332 / ((math.sqrt(1+0.5**2)*len(x_orig)) / (len(x_new)-1))

# Now use genpoly_lt.GenPoly to generate an LT file describing
# a coarse-grained DNA molecule wrapped along this curve.
# (Note: Since there is only one polymer, the "coords_multi"
#  and "name_sequence_multi" arguments contain only one list each.
#  More generally they could contain multiple lists, one for each
#  polymer in the system.)

gp = moltemplate.genpoly_lt.GenPoly()
gp.coords_multi = [x_new]

# Now specify the identity of each monomer in the polymer
# (In this case each "monomer" is a DNA base-pair, with names like "AT", "GC".)
gp.name_sequence_multi =[['AT', 'CG', 'GC', 'TA', 'AT', 'CG', 'GC', 'TA',
                          'AT', 'CG', 'GC', 'TA', 'AT', 'CG', 'GC', 'TA',
                          'AT', 'CG', 'GC', 'TA', 'AT']]

# The remaining settings are parsed from an argument list
gp.ParseArgs(['-helix', '34.2857',
              '-bond', 'Backbone', 'f', 'f',
              '-bond', 'Backbone', 'r', 'r',
              '-polymer-name', 'DNA_snippet',
              '-inherits', 'OXDNA2',
              '-header', 'import oxdna2.lt',
              '-circular', 'no'])

# Generate an .LT file and write it to the terminal (sys.stdout)
import sys
gp.WriteLTFile(sys.stdout)
