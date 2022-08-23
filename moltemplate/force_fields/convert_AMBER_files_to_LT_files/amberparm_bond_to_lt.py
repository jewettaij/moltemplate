#!/usr/bin/env python

import sys

lines_bonds = sys.stdin.readlines()
bond_style_name = 'harmonic'

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
for i in range(0, len(lines_bonds)):
    line = lines_bonds[i]
    atypes = tuple(map(str.strip, line[:6].split('-')))
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
    if (a_includes_b or (a[1] < b[1])):
        return -1
    elif (b_includes_a or (b[1] < a[1])):
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
    tokens= line.split()
    atypes = line[:6].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype1 = atype1.replace('*','star')
    atype2 = atype2.replace('*','star')
    bondtype = '@bond:'+atype1+'-'+atype2

    tokens= line[5:].split()
    keq = tokens[0]
    req = tokens[1]
    comments=' '.join(tokens[2:])
    sys.stdout.write('    bond_coeff '+bondtype+' '+bond_style_name+' '+keq+' '+req)
    if len(comments.strip()) > 0:
        sys.stdout.write('   # '+comments)
    sys.stdout.write('\n')


sys.stdout.write('  } # (end of bond_coeffs)\n')
sys.stdout.write('\n')
sys.stdout.write('  write_once("Data Bonds By Type") {\n')

for i in range(0, len(lines_sorted)):
    line = lines_sorted[i]
    atypes = line[:6].split('-')
    atype1 = atype1.replace('*','star')
    atype2 = atype2.replace('*','star')
    bondtype = '@bond:'+atype1+'-'+atype2
    #tokens= line[5:].split()
    #keq = tokens[0]
    #req = tokens[1]
    #comments=' '.join(tokens[2:])
    at1 = atype1
    at2 = atype2
    if at1 == 'X':
        at1 = '*'
    if at2 == 'X':
        at2 = '*'
    sys.stdout.write('    '+bondtype+' @atom:'+at1+' @atom:'+at2+'\n')

sys.stdout.write('  } # (end of Bonds By Type)\n')
sys.stdout.write('\n')
