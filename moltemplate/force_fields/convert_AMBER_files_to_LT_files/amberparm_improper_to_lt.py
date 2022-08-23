#!/usr/bin/env python

import sys

lines_impropers = sys.stdin.readlines()
improper_style_name = 'cvff'

# Preprocessing step.
# We want entries with more wildcards to appear earlier in the LT file.
# (Why? Moltemplate processes the rules in the "By Type" section of the
#  resulting LT file in the order in which they appear.  Putting the
#  wildcards first enables moltemplate to override generic interactions
#  containing wildcards which appear earlier in the LT file whenever it
#  finds more specific interactions that appear later in the LT file.)
# So we sort the entries ("lines") in this table so that
#  1) A line containing a string with wildcards ("X", for example "X-O2-CO-O2")
#     comes before an equivalent string without wildcards ("CO-O2-CO-O2")
#  2) Otherwise, the order is determined by their original position in the list.
lines_sorted = []
sort_keys = []
for i in range(0, len(lines_impropers)):
    line = lines_impropers[i]
    atypes = tuple(map(str.strip, line[:11].split('-')))
    sort_keys.append((atypes, i, line))
def compare_lines(a, b):
    assert((len(a) == 3) and (len(b) == 3) and (len(a[0]) == len(b[0])))
    a_includes_b = ('X' in a[0])  # does a[0] have any wildcards?
    for i in range(0, len(a[0])):
        if not (a[0][i] == 'X' or a[0][i] == b[0][i]):
            a_includes_b = False
            break
    b_includes_a = ('X' in b[0])  # does b[0] have any wildcards?
    for i in range(0, len(b[0])):
        if not (b[0][i] == 'X' or a[0][i] == b[0][i]):
            b_includes_a = False
            break
    if ((a_includes_b and (a[0] != b[0])) or (a[1] < b[1])):
        return -1
    elif ((b_includes_a and (b[0] != a[0])) or (b[1] < a[1])):
        return 1
    else:
        return 0
import functools
for _, _, line in sorted(sort_keys,
                         key = functools.cmp_to_key(compare_lines)):
    lines_sorted.append(line)


# Now loop through the list again to generate the text of the
# LT file which contains the force field parameters
sys.stdout.write('  write_once("In Settings") {\n')

for i in range(0, len(lines_sorted)):
    line = lines_sorted[i]
    atypes = line[:11].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    atype4 = atypes[3].strip()
    atype1 = atype1.replace('*','star')
    atype2 = atype2.replace('*','star')
    atype3 = atype3.replace('*','star')
    atype4 = atype4.replace('*','star')
    # NOTE: In "gaff.dat", the central atom is the third atom
    # http://archive.ambermd.org/201307/0519.html
    # We will have to select a module file from the "nbody_alt_symmetry"
    # subdirectory which uses the same convention regarding the central atom:
    impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

    tokens= line[11:].split()
    Kn = float(tokens[0])
    dn = float(tokens[1])
    n = int(float(tokens[2]))
    comments=' '.join(tokens[3:])
    if len(comments.strip()) > 0:
        comments = '    # ' + comments

    if (dn < 0.001):
        sys.stdout.write('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' 1 '+str(n)+comments+'\n')
    elif (179.999 < abs(dn) < 180.001):
        sys.stdout.write('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' -1 '+str(n)+comments+'\n')
    else:
        sys.stderr.write('Error: Illegal Improper parameters:\n'
                         '       As of 2013-8-03, LAMMPS doens hot have an improper style\n'
                         '       which can handle impropers with gamma != 0 or 180\n')
        exit(-1)



sys.stdout.write('  } # (end of improper_coeffs)\n')
sys.stdout.write('\n')
sys.stdout.write('  write_once("Data Impropers By Type (gaff_imp.py)") {\n')

for i in range(0, len(lines_sorted)):
    line = lines_sorted[i]
    atypes = line[:11].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    atype4 = atypes[3].strip()
    atype1 = atype1.replace('*','star')
    atype2 = atype2.replace('*','star')
    atype3 = atype3.replace('*','star')
    atype4 = atype4.replace('*','star')
    impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4
    at1 = atype1
    at2 = atype2
    at3 = atype3
    at4 = atype4
    if at1 == 'X':
        at1 = '*'
    if at2 == 'X':
        at2 = '*'
    if at3 == 'X':
        at3 = '*'
    if at4 == 'X':
        at4 = '*'
    sys.stdout.write('    '+impropertype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')
    #       The improper-angle is the angle between the planes
    #       defined by at1,at2,at3, and at2,at3,at3
    #       and we list the atoms in this order.
    # NOTE: In "gaff.dat", the central atom is the third atom (at3)
    #       so we have to take this into account when matching atom order.
    # http://archive.ambermd.org/201307/0519.html


sys.stdout.write('  } # (end of Impropers By Type)\n')
sys.stdout.write('\n')

# NOTE: AMBER documentation is not clear how the improper angle is defined.
#       It's not clear if we should be using the dihedral angle between
#       planes I-J-K and J-K-L.  As of 2014-4, improper_style cvff does this.
#       Even if we create improper interactions with the angle defined between
#       the wrong planes, at least the minima should be the same
#       (0 degrees or 180 degrees).
#       So I'm not too worried we are getting this detail wrong long as
#       we generate new impropers realizing that the 3rd atom (K) is the
#       central atom (according to AMBER conventions).
#
# http://structbio.vanderbilt.edu/archives/amber-archive/2007/0408.php
#
# Currently, we only apply improper torsional angles for atoms
# in a planar conformations. Is it clear?
# Junmei
