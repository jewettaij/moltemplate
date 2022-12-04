#!/usr/bin/env python3

# Author: Andrew Jewett
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2022

"""
This module contains functions (eg. ConvertAmber2Lt()) for converting AMBER
force-field files (such as FRCMOD files and DAT files) into moltemplate
format (LT files).  It can also be run as a stand-alone program.  Eg:
    amber2lt.py --in gaff.dat --out gaff.lt --name GAFF
"""

# ----------------------------------------------------------------------
# The basic atom nomenclature and conventions of AMBER files are explained here:
#   http://ambermd.org/antechamber/gaff.pdf
# For reference, the original gaff.dat file and format documentation are here:
#   http://ambermd.org/AmberTools-get.html
#   http://ambermd.org/formats.html#parm.dat
# ----------------------------------------------------------------------
# File format notes:
# How to interpret wildcard atom types ("X ") in AMBER force-field files:
#
#http://structbio.vanderbilt.edu/archives/amber-archive/2005/3444.php
#> > In the parm99 file (for example), sometimes the wild-card is used, as it
#> > is done in the following example:
#> >
#> > X -X -C -O 10.5 180. 2. JCC,7,(1986),230
#> >
#> > The first example is the specific case while the second one is the generic
#> > case. In page # 257 of the AMBER Manual, it is talking about Dihedral
#> > Angle, and how these dihedral parameters are used to calculate the
#> > energies. I am wondering what the difference between generic and specific
#> > case is for improper torsions.
#>
#> "specific" torsions are search for first, and used if a match is found. If
#> no match is found, then a search is made to see if a "generic"(aka wild-card)
#> torsion with match.
#> ...good luck...dac


import sys
import argparse
from operator import itemgetter
from collections import defaultdict
import functools


# Global variables
g_filename = __file__.split('/')[-1]
g_module_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str = '2022-9-09'
g_version_str = '0.2.1'
g_program_name = g_filename
#sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')



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



g_usage_msg = """
Typical Usage Examples:

amber2lt.py --in gaff.dat --out gaff.lt --name GAFF

(Then include the line "import gaff.lt" at the beginning of the LT
 describing the molecule(s) that needs this frcmod.  For example
#----- "benzene.lt" file -----
import "gaff.lt"
Benzene inherits GAFF {
  ...define molecule here...
}  
"""

# Example 2:
#amber2lt.py --in DABPA_AC.frcmod --out forcefield_dabpa_ac.lt --name DABPA_AC_FF



def IsFloat(s):
    try:
        x = float(s)
        return True
    except(ValueError):
        return False


def ExtractFFTableAB(lines,
                     num_atom_types,
                     allowed_section_names = None,
                     num_skip_lines = 0,
                     blanks_before_target = -1):
    count_blanks = 0
    a = -1
    b = -1
    i = num_skip_lines
    prev_line_was_blank = False
    while i < len(lines):
        line = lines[i].strip()
        if ((lines[i].strip() == "") and (not prev_line_was_blank)):
            count_blanks += 1
            if a != -1:
                b = i
                break
            else:
                prev_line_was_blank = True
                i += 1
                continue
        if ((blanks_before_target == -1) or
            (count_blanks == blanks_before_target)):
            prev_line_was_blank = False
            # select only lines which lack the pattern "??-??-??-...-?? "
            matches_pattern = len(line) > num_atom_types*3
            if matches_pattern:
                tokens = line[num_atom_types*3:].strip().split()
                if not IsFloat(tokens[0]):
                    matches_pattern = False
                    prev_line_was_blank = True
                for j in range(0, num_atom_types):
                    if line[3*j] == ' ':
                        matches_pattern = False
                for j in range(0, num_atom_types-1):
                    if line[3*j+2] != '-':
                        matches_pattern = False
                if line[num_atom_types*3-1] != ' ':
                    matches_pattern = False
            if matches_pattern:
                if a == -1:
                    a = i
            else:
                tokens = line.strip().split()
                if allowed_section_names and len(allowed_section_names) > 0:
                    if ((len(tokens) == 1) and
                        (not tokens[0] in allowed_section_names)):
                        raise InputError('Format Error: Expected ' +
                                         str(allowed_section_names) +
                                         '. Found \'' +
                                         tokens[0]+
                                         '\'. (Split file and use --bond, --angle, --dihedral, --improper arguments instead.)')
                elif len(tokens) != 0:
                    raise InputError('Format Error: Unnexpected text: "' +
                                     line.strip() + '"')
        i += 1
    if b == -1:
        b = len(lines)
    return (a, b);





def ExtractMassTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 1, ('MASS', 'MASSES'), num_skip_lines)

def ExtractBondTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 2, ('BOND', 'BONDS'), num_skip_lines)

def ExtractAngleTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 3, ('ANGL', 'ANGLE', 'ANGLES'), num_skip_lines)

def ExtractDihedralTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 4, ('DIHE', 'DIHEDRAL', 'DIHEDRALS', 'TOR', 'TORSION', 'TORSIONS'), num_skip_lines)

def ExtractImproperTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 4, ('IMPR', 'IMPROPER', 'IMPROPERS'), num_skip_lines)

def ExtractPairTextAB(lines, num_skip_lines = 0):
    return ExtractFFTableAB(lines, 1, ('NONB', 'NONBON', 'NONBOND', 'NONBONDED', 'MOD4'), num_skip_lines)



def SortLines(lines, pattern_width):
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
    for i in range(0, len(lines)):
        line = lines[i]
        # The next line converts strings like "C -CX-C8-C8" to ('C','CX','C8','C8')
        atypes = tuple(map(str.strip, line[:pattern_width].split('-')))
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
    for _, _, line in sorted(sort_keys,
                             key = functools.cmp_to_key(compare_lines)):
        lines_sorted.append(line)
    return lines_sorted





def ConvertAtomDescr2Lt(lines_mass):
    output_lines = []
    for i in range(0, len(lines_mass)):
        line = lines_mass[i]
        tokens= line.split()
        atype = tokens[0]
        mass=tokens[1]
        # what is the next number?  (the one in tokens[2]?)
        comments=' '.join(tokens[3:])
        output_lines.append('  #  @atom:'+atype+'  '+comments+'\n')
    return output_lines




def ConvertMass2Lt(lines_mass):
    output_lines = []
    output_lines.append('  write_once(\"Data Masses\") {\n')
    for i in range(0, len(lines_mass)):
        line = lines_mass[i]
        tokens= line.split()
        atype = tokens[0]
        atype = atype.replace('*','star')
        #atype = atype.replace('X','*')
        assert(atype != 'X')
        mass=tokens[1]
        # what is the next number?  (the one in tokens[2]?)
        comments=' '.join(tokens[3:])
        output_lines.append('    @atom:'+atype+' '+mass+'\n')
    output_lines.append('  } # (end of masses)\n')
    return output_lines






def ConvertBond2Lt(lines_bond):
    output_lines = []
    lines_sorted = SortLines(lines_bond, 6)
    bond_style_name = 'harmonic'
    # Now loop through the list again to generate the text of the
    # LT file which contains the force field parameters
    output_lines.append('  write_once("In Settings") {\n')
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
        output_lines.append('    bond_coeff '+bondtype+' '+bond_style_name+' '+keq+' '+req)
        if len(comments.strip()) > 0:
            output_lines.append('   # '+comments)
        output_lines.append('\n')
    output_lines.append('  } # (end of bond_coeffs)\n')
    output_lines.append('\n')
    output_lines.append('  write_once("Data Bonds By Type") {\n')
    for i in range(0, len(lines_sorted)):
        line = lines_sorted[i]
        atypes = line[:6].split('-')
        atype1 = atypes[0].strip()
        atype2 = atypes[1].strip()
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
        output_lines.append('    '+bondtype+' @atom:'+at1+' @atom:'+at2+'\n')
    output_lines.append('  } # (end of Bonds By Type)\n')
    return output_lines





def ConvertAngle2Lt(lines_angle):
    output_lines = []
    lines_sorted = SortLines(lines_angle, 8)
    angle_style_name = 'harmonic'
# Now loop through the list again to generate the text of the
# LT file which contains the force field parameters
    output_lines.append('  write_once("In Settings") {\n')
    for i in range(0, len(lines_sorted)):
        line = lines_sorted[i]
        atypes = line[:8].split('-')
        atype1 = atypes[0].strip()
        atype2 = atypes[1].strip()
        atype3 = atypes[2].strip()
        atype1 = atype1.replace('*','star')
        atype2 = atype2.replace('*','star')
        atype3 = atype3.replace('*','star')
        angletype = '@angle:'+atype1+'-'+atype2+'-'+atype3

        tokens= line[8:].split()
        keq = tokens[0]
        req = tokens[1]
        comments=' '.join(tokens[2:])
        output_lines.append('    angle_coeff '+angletype+' '+angle_style_name+' '+keq+' '+req)
        if len(comments.strip()) > 0:
            output_lines.append('   # '+comments)
        output_lines.append('\n')
    output_lines.append('  } # (end of angle_coeffs)\n')
    output_lines.append('\n')
    output_lines.append('  write_once("Data Angles By Type") {\n')
    for i in range(0, len(lines_sorted)):
        line = lines_sorted[i]
        atypes = line[:8].split('-')
        atype1 = atypes[0].strip()
        atype2 = atypes[1].strip()
        atype3 = atypes[2].strip()
        atype1 = atype1.replace('*','star')
        atype2 = atype2.replace('*','star')
        atype3 = atype3.replace('*','star')
        #tokens= line[8:].split()
        #keq = tokens[0]
        #req = tokens[1]
        #comments=' '.join(tokens[2:])
        at1 = atype1
        at2 = atype2
        at3 = atype3
        if at1 == 'X':
            at1 = '*'
        if at2 == 'X':
            at2 = '*'
        if at3 == 'X':
            at3 = '*'
        angletype = '@angle:'+atype1+'-'+atype2+'-'+atype3
        output_lines.append('    '+angletype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+'\n')
    output_lines.append('  } # (end of Angles By Type)\n')
    return output_lines









def ConvertDihedral2Lt(lines_dihedral):
    output_lines = []
    lines_sorted = SortLines(lines_dihedral, 11)
    dihedral_style_name = 'fourier'
    in_dihedral_coeffs = []
    """
       Official File Format Documentation:
             - 6 -      ***** INPUT FOR DIHEDRAL PARAMETERS *****

                        IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN

                            FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)

             IPT, ...   The atom symbols for the atoms forming a dihedral
                        angle.  If IPT .eq. 'X ' .and. LPT .eq. 'X ' then
                        any dihedrals in the system involving the atoms "JPT"
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
                    # ...I THINK (?) (Very confusing. See documentation above..)
        dn = float(tokens[2])
        n = int(float(tokens[3]))
        comments=' '.join(tokens[4:])
        if len(comments.strip()) > 0:
            comments = '    # ' + comments
        in_dihedral_coeffs.append([dihedraltype, Kn, n, dn, comments])
    # ---- processing dihedral fourier series ----
    # ---- (negative "n" values means the Fourier series is not yet complete.)
    i_orig = 1
    i = 1
    while i < len(in_dihedral_coeffs):
        type_str = in_dihedral_coeffs[i][0]
        Kn = in_dihedral_coeffs[i][1]
        n = in_dihedral_coeffs[i][2]
        dn = in_dihedral_coeffs[i][3]
        comments = in_dihedral_coeffs[i][-1]
        #sys.stderr.write("orig_dihedral_coeffs["+str(i_orig)+"] = "+str(in_dihedral_coeffs[i])+"\n")
        if ((i>0) and (in_dihedral_coeffs[i-1][-3] < 0)):
            assert(in_dihedral_coeffs[i-1][0] == in_dihedral_coeffs[i][0])
            in_dihedral_coeffs[i-1][-3] = -in_dihedral_coeffs[i-1][-3]
            old_comments = in_dihedral_coeffs[i-1][-1]
            if (len(old_comments) > 0) and (len(comments) == 0):
                comments = old_comments
            in_dihedral_coeffs[i-1][-1] = Kn
            in_dihedral_coeffs[i-1].append(n)
            in_dihedral_coeffs[i-1].append(dn)
            in_dihedral_coeffs[i-1].append(comments)
            del in_dihedral_coeffs[i]
        else:
            #sys.stderr.write("in_dihedral_coeffs["+str(i-1)+"] = "+str(in_dihedral_coeffs[i-1])+"\n")
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
    output_lines.append('  write_once(\"In Settings\") {\n    ')
    output_lines.append('\n    '.join(in_dihedral_coeffs)+'\n')
    output_lines.append('  } # (end of dihedral_coeffs)\n')
    output_lines.append('\n')
    output_lines.append('  write_once("Data Dihedrals By Type") {\n')
    for i in range(0, len(lines_sorted)):
        line = lines_sorted[i]
        atypes = line[:11].split('-')
        if (i > 0) and (atypes == atypes_prev):
            continue   # no need to redundantly repeat the same rule multiple times
        atypes_prev = atypes
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
        output_lines.append('    '+dihedraltype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')
    output_lines.append('  } # (end of Dihedrals By Type)\n')
    return output_lines







def ConvertImproper2Lt(lines_improper):
    output_lines = []
    lines_sorted = SortLines(lines_improper, 11)
    improper_style_name = 'cvff'
    # Now loop through the list again to generate the text of the
    # LT file which contains the force field parameters
    output_lines.append('  write_once("In Settings") {\n')
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
        # subdirectory which uses the same convention regarding the central atom
        impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4
        tokens= line[11:].split()
        Kn = float(tokens[0])
        dn = float(tokens[1])
        n = int(float(tokens[2]))
        comments=' '.join(tokens[3:])
        if len(comments.strip()) > 0:
            comments = '    # ' + comments
        if (dn < 0.001):
            output_lines.append('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' 1 '+str(n)+comments+'\n')
        elif (179.999 < abs(dn) < 180.001):
            output_lines.append('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' -1 '+str(n)+comments+'\n')
        else:
            sys.stderr.write('Error: Illegal Improper parameters:\n'
                             '       As of 2013-8-03, LAMMPS doens hot have an improper style\n'
                             '       which can handle impropers with gamma != 0 or 180\n')
            exit(-1)
    output_lines.append('  } # (end of improper_coeffs)\n')
    output_lines.append('\n')
    output_lines.append('  write_once("Data Impropers By Type (gaff_imp.py)") {\n')
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
        output_lines.append('    '+impropertype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')
        #       The improper-angle is the angle between the planes
        #       defined by at1,at2,at3, and at2,at3,at3
        #       and we list the atoms in this order.
        # NOTE: In "gaff.dat", the central atom is the third atom (at3)
        #       so we have to take this into account when matching atom order.
        # http://archive.ambermd.org/201307/0519.html
    output_lines.append('  } # (end of Impropers By Type)\n')
    return output_lines
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








def ConvertPair2Lt(lines_pair,
                   lines_mass=[]):
    output_lines = []
    pair_style = 'lj/charmm/coul/long'
    atom_types_in_pairs = set([])
    # NOTE: LAMMPS complains if you attempt to use lj/charmm/coul/long
    # on a system if it does not contain any charged particles.
    # Moltemplate does not assign atomic charge,
    # so this problem occurs frequently unless the user remembers to
    # supply the charge.
    output_lines.append('  write_once(\"In Settings\") {\n')
    for i in range(0, len(lines_pair)):
        line = lines_pair[i]
        tokens= line.split()
        atype = tokens[0]
        atype = atype.replace('*','star')
        #atype = atype.replace('X','*')
        assert(atype != 'X')
        atom_types_in_pairs.add(atype)
        # UGGHHH
        # OLD CODE:
        #sig=tokens[1]
        #   CORRECTION #1
        # It looks the number in this part of the file is an atom radii, not a
        # diameter.  In other words, this number is 0.5*sigma instead of sigma.
        # So we multiply it by 2.0.
        #sig=str(2.0*float(tokens[1]))
        #
        #   CORRECTION #2
        # It also appears as though they are using this convention for LennardJones
        # U(r)=epsilon*((s/r)^12-2*(s/r)^6)    instead of  4*eps*((s/r)^12-(s/r)^6)
        #    ...where "s" is shorthand for "sigma"..
        # This means we must ALSO multiply sigma in gaff.dat by 2**(-1.0/6)
        # (This change makes the two U(r) formulas equivalent.)
        #  I had to figure this out by iterations of trial and error.
        #  The official AMBER documentation is quite vague about the LJ parameters.
        #  -Andrew 2014-5-19
        # http://ambermd.org/formats.html#parm.dat
        # http://structbio.vanderbilt.edu/archives/amber-archive/2009/5072.php)
        sig=str(float(tokens[1])*2.0*pow(2.0, (-1.0/6.0)))
        eps=tokens[2]
        comments=' '.join(tokens[3:])
        output_lines.append('    pair_coeff @atom:'+atype+' @atom:'+atype+' '+pair_style+' '+eps+' '+sig)
        if len(comments.strip()) > 0:
            output_lines.append('   # '+comments)
        output_lines.append('\n')
    # Any missing pair_coeffs (for atoms that are defined in the Masses section)
    # should be supplied with default values (so that LAMMPS doesn't complain).
    for i in range(0, len(lines_mass)):
        line = lines_mass[i]
        tokens= line.split()
        atype = tokens[0]
        if not (atype in atom_types_in_pairs):
            eps = '0.0'
            sig = '1.0'
            output_lines.append('    pair_coeff @atom:'+atype+' @atom:'+atype+' '+pair_style+' '+eps+' '+sig+' # (default parameters)\n')
            sys.stderr.write('WARNING: No Lennard-Jones parameters defined for atom type "'+atype+'" in source file.\n')
    output_lines.append('  } # (end of pair_coeffs)\n')
    return output_lines







def ConvertAmberSections2Lt(object_name,
                            lines_mass,
                            lines_bond,
                            lines_angle,
                            lines_dihedral,
                            lines_improper,
                            lines_pair,
                            lines_preamble = []):
    lines_atom_descr_new = ConvertAtomDescr2Lt(lines_mass)
    lines_mass_new       = ConvertMass2Lt(lines_mass)
    lines_bond_new       = ConvertBond2Lt(lines_bond)
    lines_angle_new      = ConvertAngle2Lt(lines_angle)
    lines_dihedral_new   = ConvertDihedral2Lt(lines_dihedral)
    lines_improper_new   = ConvertImproper2Lt(lines_improper)
    lines_pair_new       = ConvertPair2Lt(lines_pair, lines_mass)

    lines_prefix = []
    if object_name and object_name != "":
        lines_prefix.append(object_name+' {\n')
    lines_suffix = []
    if object_name and object_name != "":
        lines_suffix.append('} # '+object_name+'\n')

    init_section = """
  write_once("In Init") {
    # Default styles and settings for AMBER based force-fields:
    units           real
    atom_style      full
    bond_style      hybrid harmonic
    angle_style     hybrid harmonic
    dihedral_style  hybrid fourier
    improper_style  hybrid cvff
    pair_style      hybrid lj/charmm/coul/long 9.0 10.0 10.0
    kspace_style    pppm 0.0001

    # NOTE: If you do not want to use long-range coulombic forces,
    #       comment out the two lines above and uncomment this line:
    # pair_style      hybrid lj/charmm/coul/charmm 9.0 10.0

    pair_modify     mix arithmetic
    special_bonds   amber
  }
"""

    lines_init = [init_section]

    if lines_preamble and len(lines_preamble) > 0:
        lines_preamble += ['\n']
    else:
        lines_preamble = []

    output_lines = (lines_preamble +
                    lines_prefix + ['\n'] +
                    lines_atom_descr_new + ['\n'] +
                    lines_mass_new + ['\n'] +
                    lines_pair_new + ['\n'] +
                    lines_bond_new + ['\n'] +
                    lines_angle_new + ['\n'] +
                    lines_dihedral_new + ['\n'] +
                    lines_improper_new + ['\n'] +
                    lines_init + ['\n'] +
                    lines_suffix + ['\n'])

    return output_lines




def ConvertAmber2Lt(object_name, lines):
    ip = 1
    a, b = ExtractMassTextAB(lines, ip)
    lines_mass = lines[a:b]
    ip = b
    a, b = ExtractBondTextAB(lines, ip)
    lines_bond = lines[a:b]
    ip = b
    a, b = ExtractAngleTextAB(lines, ip)
    lines_angle = lines[a:b]
    ip = b
    a, b = ExtractDihedralTextAB(lines, ip)
    lines_dihedral = lines[a:b]
    ip = b
    a, b = ExtractImproperTextAB(lines, ip)
    lines_improper = lines[a:b]
    ip = b
    a, b = ExtractPairTextAB(lines, ip)
    lines_pair = lines[a:b]
    ip = b
    return ConvertAmberSections2Lt(object_name,
                                   lines_mass,
                                   lines_bond,
                                   lines_angle,
                                   lines_dihedral,
                                   lines_improper,
                                   lines_pair)






def main():
    # Inform the user what version of the software they are using
    sys.stderr.write(g_program_name + ' v' +
                     g_version_str + ' ' + g_date_str + '\n')
    sys.stderr.write('WARNING: THIS IS EXPERIMENTAL SOFTWARE (2022-8-24)\n')
    try:
        # Now parse the arguments passed to the program (if any)
        ap = argparse.ArgumentParser()
        ap.add_argument('-i', '-in', '--in',
                        dest='in_filename',
                        required=False,
                        help='name of the .FRCMOD or .DAT file you want to convert (eg "gaff2.lt". if unspecified, stdin is used)')
        ap.add_argument('-o', '-out', '--out',
                        dest='out_filename',
                        required=False,
                        help='name of the LT file you want to create (if unspecified, stdout is used)')
        ap.add_argument('-name', '--name',
                        dest='object_name',
                        required=False,
                        help='name of the object (force field) you want to create (eg. "GAFF2")')
        ap.add_argument('-mass', '-masses','--mass', '--masses',
                        dest='filename_mass',
                        required=False,
                        help='name of a file containing an excerpt with only mass information')
        ap.add_argument('-bond', '-bonds','--bond', '--bonds',
                        dest='filename_bond',
                        required=False,
                        help='name of a file containing an excerpt with only bond information')
        ap.add_argument('-angle', '-angles','--angle', '--angles',
                        dest='filename_angle',
                        required=False,
                        help='name of a file containing an excerpt with only (3-body) angle information')
        ap.add_argument('-dihedral', '-dihedrals','--dihedral', '--dihedrals',
                        dest='filename_dihedral',
                        required=False,
                        help='name of a file containing an excerpt with only (4-body) dihedral information')
        ap.add_argument('-improper', '-impropers','--improper', '--impropers',
                        dest='filename_improper',
                        required=False,
                        help='name of a file containing an excerpt with only (4-body) improper information')
        ap.add_argument('-pair', '-pairs','--pair', '--pairs',
                        dest='filename_pair',
                        required=False,
                        help='name of a file containing an excerpt with only pair (non-bonded) information')

        args = ap.parse_args()

        # Figure out the names of the LT file the user wants to create
        # (By default, this program will write to the terminal (sys.stdout).)
        if args.out_filename:
            out_file = open(args.out_filename, 'w')
        else:
            out_file = sys.stdout

        # Now figure out the names of the file(s) the user wants us to read

        if ((args.filename_mass and args.filename_mass != '') and
            (args.filename_bond and args.filename_bond != '') and
            (args.filename_angle and args.filename_angle != '') and
            (args.filename_dihedral and args.filename_dihedral != '') and
            (args.filename_improper and args.filename_improper != '') and
            (args.filename_pair and args.filename_pair != '')):
            with open(finename_mass, 'r') as f:
                lines_mass = f.readlines()
            with open(finename_bond, 'r') as f:
                lines_bond = f.readlines()
            with open(finename_angle, 'r') as f:
                lines_angle = f.readlines()
            with open(finename_dihedral, 'r') as f:
                lines_dihedral = f.readlines()
            with open(finename_improper, 'r') as f:
                lines_improper = f.readlines()
            with open(finename_pair, 'r') as f:
                lines_pair = f.readlines()
        else:
            # (By default, this program will read from the terminal (sys.stdin))
            if args.in_filename:
                in_file = open(args.in_filename, 'r')
            else:
                in_file = sys.stdin
            lines = in_file.readlines()
            ip = 1  # keep track of how many lines we have parsed so far
            if args.in_filename:
                in_file.close()
            # Extract the mass information
            a, b = ExtractMassTextAB(lines, ip)
            lines_mass = lines[a:b]
            ip = b
            # Extract the bond information
            a, b = ExtractBondTextAB(lines, ip)
            lines_bond = lines[a:b]
            ip = b
            # Extract the angle information
            a, b = ExtractAngleTextAB(lines, ip)
            lines_angle = lines[a:b]
            ip = b
            # Extract the dihedral information
            a, b = ExtractDihedralTextAB(lines, ip)
            lines_dihedral = lines[a:b]
            ip = b
            # Extract the improper information
            a, b = ExtractImproperTextAB(lines, ip)
            lines_improper = lines[a:b]
            ip = b
            # Extract the pair information
            a, b = ExtractPairTextAB(lines, ip)
            lines_pair = lines[a:b]
            ip = b

        preamble = \
"""####################################################################
# To use this, LAMMPS must be compiled with the EXTRA-MOLECULE package
# (See here for details: https://docs.lammps.org/Build_package.html)
#######################################################################
#    This moltemplate (LT) file was generated automatically using
""" \
    + '#      amber2lt.py ' + ' '.join(sys.argv[1:]) + \
"""
#######################################################################
#      USAGE SUGGESTIONS:
# Suppose you named your force field "GAFF"  (using amber2lt.py's "--name"
# argument) and saved it in a file named "gaff.lt".  Then, in order to create
# molecules that use this force field, create a file for each type of molecule
# (eg "benzene.lt") and use this format:
# #---- "benzene.lt" file -----
# import "gaff.lt"
# Benzene inherits GAFF {
#   write('Data Atoms') {
#     $atom:c1  $mol  @atom:ca  -0.115   -0.739   1.189  -0.00733
#     $atom:c2  $mol  @atom:ca  -0.115    0.614   1.208   0.35167
#        :        :        :
#   }
#   write("Data Bond List") {
#     $bond:b1 $atom:C1 $atom:C2
#        :
#   }
# }
#
# ------- How to generate molecule files (such as "benzene.lt") -------
#
# You can try to generate these files manually, but you must be careful
# to choose the correct @atom types for each atom (eg "@atom:ca"),
# and you must obtain the atomic charge (column 4) by other means. (See below.)
#
# Recommended method:
#
# Use AmberTools to create a MOL2 file containing one of the molecules
# you want to simulate.  Suppose you use AmberTools to prepare a MOL2 file
# containing a single benzene molecule (named "benzene.mol2")
# and a corresponding FRCMOD file (named "benzene.frcmod").  Use:
#
# mol22lt.py --in benzene.mol2 --out benzene.lt --name Benzene \\
#            --ff MyForceField --ff-file my_force_field.lt
# amber2lt.py --in benzene.frcmod --name MyForceField >> my_force_field.lt
#
# ...The resulting "benzene.lt" file will have this format:
# #---- "benzene.lt" file -----
# import "my_force_field.lt"
# Benzene inherits MyForceField {
#   :
# }
# Repeat this for each different type of molecule you want in your simulation.
# (If it is a molecule with multiple chains, do this for every chain.)
# (If you want to start over, delete the "my_force_field.lt" file.)
#######################################################################

"""
        output_lines = ConvertAmberSections2Lt(args.object_name,
                                               lines_mass,
                                               lines_bond,
                                               lines_angle,
                                               lines_dihedral,
                                               lines_improper,
                                               lines_pair,
                                               [preamble])

        out_file.write(''.join(output_lines))


    except (ValueError, InputError) as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(-1)


if __name__ == '__main__':
    main()
