#!/usr/bin/env python

# SOME UGLY CODE HERE

import sys

lines_dihedrals = sys.stdin.readlines()
dihedral_style_name = 'fourier'
in_dihedral_coeffs = []

# Preprocessing step.
# We want entries with more wildcards to appear earlier in the LT file.
# (Why? Moltemplate processes the rules in the "By Type" section of the
#  resulting LT file in the order in which they appear.  Putting the
#  wildcards first enables moltemplate to override generic interactions
#  containing wildcards which appear earlier in the LT file whenever it
#  finds more specific interactions that appear later in the LT file.)
# So we sort the entries ("lines") in this table according to
#  1) the number of wildcard atom types (ie. atom types named "X")
#  2) the original order the atoms appear in the file
#     so that we can preserve the original order whenever possible
lines_sorted = []
sort_keys = []
for i in range(0, len(lines_dihedrals)):
    line = lines_dihedrals[i]
    # The next line converts strings like "C -CX-C8-C8" to ('C','CX','C8','C8')
    atypes = tuple(map(str.strip, line[:11].split('-')))
    sort_keys.append((atypes, i, line))
def compare_lines(a, b):
    assert((len(a) == 3) and (len(b) == 3) and (len(a[0]) == len(b[0])))
    if (a[0] == b[0]):
        if a[1] < b[1]:
            return -1
        elif a[1] > b[1]:
            return 1
        else:
            return 0
    a_includes_b = True
    for i in range(0, len(a[0])):
        if not (a[0][i] == 'X' or a[0][i] == b[0][i]):
            a_includes_b = False
            break
    b_includes_a = True
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
for i in range(0, len(lines_sorted)):
    line   = lines_sorted[i]
    atypes = line[:11].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    atype4 = atypes[3].strip()
    atype1 = atype1.replace('*','star')
    atype2 = atype2.replace('*','star')
    atype3 = atype3.replace('*','star')
    atype4 = atype4.replace('*','star')
    dihedraltype = '@dihedral:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

    tokens= line[11:].split()
    npth = float(tokens[0])
    Kn = float(tokens[1])
    Kn /= npth  # The coeff for each fourier term is Kn/npth
                # ...I THINK (?).  (Very confusing.  See documentation below...)
    dn = float(tokens[2])
    n = int(float(tokens[3]))
    comments=' '.join(tokens[4:])
    #print("dihedraltype = "+str(dihedraltype)+
    #      ", tokens["+str(i)+"] = "+str(tokens))
    #print("n = "+str(n))
    if len(comments.strip()) > 0:
        comments = '    # ' + comments
    in_dihedral_coeffs.append([dihedraltype, Kn, n, dn, comments])
    #print(Kn, n, dn)


#for entry in in_dihedral_coeffs:
#    print(entry)
#exit()


# ---- processing dihedral fourier series ----
# ---- (negative "n" values means the Fourier series is not yet complete.)

i_orig = 0
i = 1
while i < len(in_dihedral_coeffs):
    type_str = in_dihedral_coeffs[i][0]
    Kn = in_dihedral_coeffs[i][1]
    n = in_dihedral_coeffs[i][2]
    dn = in_dihedral_coeffs[i][3]
    comments = in_dihedral_coeffs[i][-1]

    #print("orig_dihedral_coeffs["+str(i_orig)+"] = "+str(in_dihedral_coeffs[i]))
    #if (i>0):
    #    sys.stderr.write('prev_n='+str(in_dihedral_coeffs[i-1][-3])+'\n')
    #sys.stderr.write('n='+str(n)+'\n')

    if ((i>0) and (in_dihedral_coeffs[i-1][-3] < 0)):

        #sys.stdout.write('interaction_before_append: '+str(in_dihedral_coeffs[i-1])+'\n')
        #if not (in_dihedral_coeffs[i-1][0] == in_dihedral_coeffs[i][0]):
        #    print(''.join(lines_sorted))
        #    print(' '.join(map(str,in_dihedral_coeffs[i-1])))
        #    print(' '.join(map(str,in_dihedral_coeffs[i])))
        #    print("i="+str(i)+
        #          ", in_dihedral_coeffs[i-1][0] == "+str(in_dihedral_coeffs[i-1][0])+
        #          ", in_dihedral_coeffs[i][0] == "+str(in_dihedral_coeffs[i][0]))
        #    print("in_dihedral_coeffs[i-1][2]="+str(in_dihedral_coeffs[i-1][2]))
        #    print("in_dihedral_coeffs[i][2]="+str(in_dihedral_coeffs[i][2]))
        assert(in_dihedral_coeffs[i-1][0] == in_dihedral_coeffs[i][0])
        in_dihedral_coeffs[i-1][-3] = -in_dihedral_coeffs[i-1][-3]
        old_comments = in_dihedral_coeffs[i-1][-1]
        if (len(old_comments) > 0) and (len(comments) == 0):
            comments = old_comments
        in_dihedral_coeffs[i-1][-1] = Kn
        in_dihedral_coeffs[i-1].append(n)
        in_dihedral_coeffs[i-1].append(dn)
        in_dihedral_coeffs[i-1].append(comments)
        #sys.stdout.write('interaction_after_append: '+str(in_dihedral_coeffs[i-1])+'\n')
        del in_dihedral_coeffs[i]
    else:
        #print("in_dihedral_coeffs["+str(i-1)+"] = "+str(in_dihedral_coeffs[i-1]))
        i += 1
    i_orig += 1


for i in range(0, len(in_dihedral_coeffs)):
    type_str = in_dihedral_coeffs[i][0]
    params = in_dihedral_coeffs[i][1:]
    params = list(map(str, params))
    num_fourier_terms = int((len(params)-1)/3)
    dihedral_coeff_str = 'dihedral_coeff '+type_str+' '+\
        dihedral_style_name+' '+ \
        str(num_fourier_terms)+' '+ \
        ' '.join(params)
    in_dihedral_coeffs[i] = dihedral_coeff_str

# ---- finished processing dihedral fourier series ----


sys.stdout.write('  write_once(\"In Settings\") {\n    ')
sys.stdout.write('\n    '.join(in_dihedral_coeffs)+'\n')
sys.stdout.write('  } # (end of dihedral_coeffs)\n')





sys.stdout.write('\n')

sys.stdout.write('  write_once("Data Dihedrals By Type") {\n')

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
    dihedraltype = '@dihedral:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4
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
    sys.stdout.write('    '+dihedraltype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')

sys.stdout.write('  } # (end of Dihedrals By Type)\n')
sys.stdout.write('\n')


"""
         - 6 -      ***** INPUT FOR DIHEDRAL PARAMETERS *****

                    IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN

                        FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)

         IPT, ...   The atom symbols for the atoms forming a dihedral
                    angle.  If IPT .eq. 'X ' .and. LPT .eq. 'X ' then
                    any dihedrals in the system involving the atoms "JPT" and
                    and "KPT" are assigned the same parameters.  This is
                    called the general dihedral type and is of the form
                    "X "-"JPT"-"KPT"-"X ".

         IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

         PK         The barrier height divided by a factor of 2.

         PHASE      The phase shift angle in the torsional function.

                    The unit is degrees.

         PN         The periodicity of the torsional barrier.
                    NOTE: If PN .lt. 0.0 then the torsional potential
                          is assumed to have more than one term, and the
                          values of the rest of the terms are read from the
                          next cards until a positive PN is encountered.  The
                          negative value of pn is used only for identifying
                          the existence of the next term and only the
                          absolute value of PN is kept.

            The input is terminated by a blank card.
"""
