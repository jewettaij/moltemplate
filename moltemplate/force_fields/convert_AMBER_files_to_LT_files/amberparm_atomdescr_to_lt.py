#!/usr/bin/env python

import sys

lines_gaff = sys.stdin.readlines()

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    tokens= line.split()
    atype = tokens[0]
    mass=tokens[1]
    # what is the next number?  (the one in tokens[2]?)
    comments=' '.join(tokens[3:])
    sys.stdout.write('#  @atom:'+atype+'  '+comments+'\n')

sys.stdout.write('\n')
