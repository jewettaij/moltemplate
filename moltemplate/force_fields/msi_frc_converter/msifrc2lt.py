#! /usr/bin/env python

"""
This standalone python script can be used to convert the force-fields in MSI
format (.FRC files, a.k.a. "BIOSYM", "DISCOVERY" format)
...into MOLTEMPLATE/LAMMPS format (.LT format).

Once converted into moltemplate (.LT) format, users can use these files with 
MOLTEMPLATE to prepare LAMMPS simulations of molecules using these force fields 
(without needing any additional software such as msi2lmp).

There are several examples of MSI files in the "tools/msi2lmp/frc_files/"
directory which is distributed with LAMMPS.

Limitations:

Currently (2017-2) this script ignores the "template" information in .FRC files.
When defining a new type of molecule, the user must carefully choose the
complete atom type for each type of atom in the molecule.  In other words,
MOLTEMPLATE will not attempt to determine (from local context) whether 
a carbon atom somewhere in your molecule happens to be an SP3 carbon 
(ie. "c4" in the COMPASS force-field), or an aromatic carbon ("c3a"), 
or something else (for example).  This information is typically contained 
in the "templates" section of these files, and this script currently ignores
that information.  Instead, the user must determine which type of carbon atom
it is manually, for all of the carbon atoms in that kind of molecule.
(This only needs to be done once per molecule definition.
 Once a type of molecule is defined, it can be copied indefinitely.)

"""


__author__ = 'Andrew Jewett'
__version__ = '0.1.2'
__date__ = '2017-3-16'


import sys
import os

from collections import defaultdict, OrderedDict
from operator import itemgetter
from math import *

g_program_name = __file__.split('/')[-1]


doc_msg = \
    "Typical Usage:\n\n" + \
    "   " + g_program_name + " -name COMPASS < compass_published.frc > compass.lt\n\n" + \
    "   where \"compass_published.frc\" is a force-field file in MSI format.\n" + \
    "         \"comass.lt\" is the corresponding file converted to moltemplate format\n" + \
    "   and   \"COMPASS\" is the name that future moltemplate users will use to refer\n" + \
    "         to this force-field (optional).\n" + \
    "Optional Arguments\n" + \
    "   -name FORCEFIELDNAME # Give the force-field a name\n" + \
    "   -file FILE_NAME      # Read force field parameters from a file\n" + \
    "   -url URL             # Read force field parameters from a file on the web\n" + \
    "   -atoms \"QUOTED LIST\" # Restrict output to a subset of atom types\n" + \
    "  Sometimes an FRC file contains multiple versions.  In that case,\n"+\
    "  you can select between them using these optional arguments:\n"+\
    "   -pair-style \"PAIRSTYLE ARGS\" # LAMMPS pair style and cutoff arg(s)\n" + \
    "   -bond-style BONDSTYLE  # desired LAMMPS bond style (default: \"class2\")\n" + \
    "   -angle-style ANGLESTYLE  # desired LAMMPS angle style\n" + \
    "   -dihedral-style DIHEDRALSTYLE  # desired LAMMPS dihedral style\n" + \
    "   -improper-style IMPROPERSTYLE  # desired LAMMPS improper style\n" + \
    "   -hbond-style \"HBONDTYLE ARGS\" # LAMMPS hydrogen-bond style and args\n"


#   "   -auto                # Consider auto_equivalences in the .frc file \n"+\



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


def NSplitQuotedString(string,
                       nmax,
                       quotes,
                       delimiters=' \t\r\f\n',
                       escape='\\',
                       comment_char='#'):
    """
    Split a quoted & commented string into at most "nmax" tokens (if nmax>0),
    where each token is separated by one or more delimeter characters
    in the origingal string, and quoted substrings are not split,
    This function returns a list of strings.  Once the string is split Nmax
    times, any remaining text will be appended to the last entry of the list.
    Comments are stripped from the string before splitting begins.
    """
    tokens = []
    token = ''
    reading_token = True
    escaped_state = False
    quote_state = None
    for c in string:

        if (c in comment_char) and (not escaped_state) and (quote_state == None):
            if len(token) > 0:
                tokens.append(token)
            return tokens

        elif (c in delimiters) and (not escaped_state) and (quote_state == None):
            if reading_token:
                if (nmax == 0) or (len(tokens) < nmax-1):
                    if len(token) > 0:
                        tokens.append(token)
                    token = ''
                    reading_token = False
                else:
                    token += c
        elif c in escape:
            if escaped_state:
                token += c
                reading_token = True
                escaped_state = False
            else:
                escaped_state = True
                # and leave c (the '\' character) out of token
        elif (c in quotes) and (not escaped_state):
            if (quote_state != None):
                if (c == quote_state):
                    quote_state = None
            else:
                quote_state = c
            token += c
            reading_token = True
        else:
            if (c == 'n') and (escaped_state == True):
                c = '\n'
            elif (c == 't') and (escaped_state == True):
                c = '\t'
            elif (c == 'r') and (escaped_state == True):
                c = '\r'
            elif (c == 'f') and (escaped_state == True):
                c = '\f'
            token += c
            reading_token = True
            escaped_state = False

    if len(token) > 0:
        tokens.append(token)
    return tokens




def SplitQuotedString(string,
                      quotes='\'\"',
                      delimiters=' \t\r\f\n',
                      escape='\\',
                      comment_char='#'):

    return NSplitQuotedString(string,
                              0,
                              quotes,
                              delimiters,
                              escape,
                              comment_char)




def RemoveOuterQuotes(text, quotes='\"\''):
    if ((len(text) >= 2) and (text[0] in quotes) and (text[-1] == text[0])):
        return text[1:-1]
    else:
        return text


def ReverseIfEnds(l_orig):
    """
    Convenient to have a one-line macro for swapping list order if first>last
    """
    l = [x for x in l_orig]
    if l[0] > l[-1]:
        l.reverse()
    return l



#def Repl(tokens, a, b):
#    return [(b if x==a else x) for x in tokens]


def EncodeAName(s):
    """
    Handle * characters in MSI atom names
    """

    # If the atom name begins with *, then it is a wildcard
    if s[:1] == '*': # special case: deal with strings like  *7
        return 'X'   # These have special meaning.  Throw away the integer.
                     # (and replace the * with an X)

    # If the * character occurs later on in the atom name, then it is actually
    # part of the atom's name.  (MSI force fields use many strange characters in
    # atom names.)  Here we change the * to \* to prevent the atom name from
    # being interpreted as a wild card in the rules for generating bonds,
    # angles, dihedrals, and impropers.
    
    return s.replace('*','\\*')  # this prevents ttree_lex.MatchesAll()
                                 # from interpreting the '*' as a wildcard

    # alternately:
    #return s.replace('*','star') # '*' is reserved for wildcards in moltemplate
    #                             # 'star' is a string that is unused in any 
    #                             # of the force fields I have seen so far.
    

def DetermineAutoPriority(anames):
    #scan through list of strings anames, looking for patterns of the form
    #*n
    #where n is an integer.
    #(These patterns are used by MSI software when using "auto_equivalences"
    # to look up force field parameters for bonded interactions.)
    #Make sure this pattern only appears once and return n to the caller.
    n = None
    for i in range(0, len(anames)):
        if anames[:1] == '*':
            if n == None:
                n = int(anames[1:])
            elif n != int(anames[1:]):
                raise InputError('Error: Inconsistent priority integers in the following interaction:\n'
                                 '      ' + ' '.join(anames) + '\n')
    if n == None:
        return 0
    else:
        return n


def DeterminePriority(is_auto,
                      anames,
                      version):
    """
    Determine the priority of an interaction from 
    1) whether or not it is an "auto" interaction
    2) what is the force-field "version" (a number)
    3) what are the names of the atoms (for auto_equivalences only,
       some atom "names" are wildcards followed by integers. use the integer)
    """

    if is_auto:
        n = DetermineAutoPriority(anames)
        return (is_auto, n)
    else:
        return (is_auto, version)

def IsAutoInteraction(interaction_name):
    return interaction_name.find('auto') == 0

def EncodeInteractionName(anames,
                          is_auto = None):
    if is_auto == None:
        is_auto = False
        # Is the line containing anames from an "_auto" section of 
        # the FRC file?  (I am trying to infer this from whether or 
        # not any of the atom names are followed by the '_' character.)
        for s in anames:  
            if ((len(s)>0) and (s[-1] == '_')):
                is_auto = True
    if is_auto:
        priority = DetermineAutoPriority(anames)
        # (If an atom name is a wildcard '*' followed by 
        #  an integer, DetermineAutoPriority() will return 
        #  that integer.  Otherwise it will return '')
        return 'auto' + str(priority)+','.join(anames)
    return ','.join(anames)

def ExtractANames(interaction_name):
    if IsAutoInteraction(interaction_name):
        return interaction_name[5:].split(',')
    return interaction_name.split(',')


def OOPImproperNameSort(aorig):
    assert(len(aorig) == 4)
    atom_names = map(EncodeAName, aorig)
    if atom_names[0] < atom_names[3]:
        return (atom_names, [0,1,2,3])
    else:
        return ([atom_names[3],
                 atom_names[1],
                 atom_names[2],
                 atom_names[0]],
                [3,1,2,0])


def Class2ImproperNameSort(aorig):
    """
    This function takes a list of 4 strings as an argument representing 4 atom
    names for atoms participating in an "improper" ("wilson-out-of-plane")
    interaction.  This function assumes the second atom is the central ("hub") 
    atom in the interaction, and it sorts the remaining atoms names.
    This function also replaces any occurence of \"*\" with \"X\".
    The new list is returned to the caller, along with the permutation.
    """
    assert(len(aorig) == 4)
    atom_names = map(EncodeAName, aorig)
    z = zip([atom_names[0], atom_names[2], atom_names[3]],
            [0,2,3])
    z.sort()
    l = [z[0][0], atom_names[1], z[1][0], z[2][0]]
    p = [z[0][1], 1, z[1][1], z[2][1]]
    return (l, p)


def ImCrossTermID(atom_names):
    if atom_names[0] <= atom_names[3]:
        cross_name = (atom_names[0]+','+atom_names[1]+','+
                      atom_names[2]+','+atom_names[3])





def DoAtomsMatchPattern(anames, pattern):
    """ 
    Check whether the list of atom names "anames" matches "pattern"
    (Both arguments are lists of strings, but some of the strings 
    in pattern may contain wildcard characters followed by 
    "priority" numbers.  Matches with lower priority numbers are
    given preference whenever multiple distinct matches are found.
    (Note: This function does not check patterns in reverse order.)
    """
    #sys.stderr.write('DEBUG: checking whether '+str(anames)+' matches '+str(pattern)+'\n')
    assert(len(anames) == len(pattern))
    matched = True
    for d in range(0, len(pattern)):
        if (pattern[d] == anames[d]) or (pattern[d][0] == '*'):
            if pattern[d][0] == '*':
                priority = int(pattern[d][1:])
            else:
                priority = 0
        else:
            matched = False
    if matched:
        #sys.stderr.write('DEBUG: '+str(anames)+' matches '+str(pattern)+'\n')
        return priority
    else:
        return None


def LookupBondLength(a1, a2,
                     atom2equiv_bond,
                     bond2r0,
                     atom2auto_bond,
                     bond2r0_auto):
    """ 
    Try to find bond parameters between atoms whose original
    atom names (without equivalences) are a1 and a2.
    Then return both the equilibrium bond length for that bond,
    as well as the equivalent atom names used to lookup that bond.
    (These could be stored in either atom2equiv_bond or atom2auto_bond.)
    If a match was not found, return None.
    """
    return_val = None
    anames = (atom2equiv_bond[a1], atom2equiv_bond[a2])
    bond_name = EncodeInteractionName(ReverseIfEnds(anames))
    if bond_name in bond2r0:
        return_val = (bond2r0[bond_name], [anames[0], anames[1]])
    # If no bond between these atoms is defined, 
    # check the bonds in the _auto section(s)
    # This is a lot messier.
    elif ((a1 in atom2auto_bond) and (a2 in atom2auto_bond)):
        anames = [atom2auto_bond[a1], atom2auto_bond[a2]]
        # Because _auto interactions can contain wildcards,
        # there can be multiple entries in bond2r0_auto[]
        # for the same list of atom names, and we have to
        # consider all of them, and pick the one with the
        # most priority (ie. whose priority number is lowest).
        # (Note: The MSI file format uses low priority numbers
        #  to indicate high priority.  Somewhat confusing.)
        HUGE_VAL = 2000000000
        best_priority = HUGE_VAL
        pattern = ['','']
        for (pattern[0],pattern[1]), r0 in bond2r0_auto.items():
            priority = DoAtomsMatchPattern(anames, pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (r0, [anames[0], anames[1]])
            anames.reverse() # now check of the atoms in reverse order match
            priority = DoAtomsMatchPattern(anames, pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (r0, [anames[1], anames[0]]) #preserve atom order
        #if return_val != None:
        #    sys.stderr.write('DEBUG: For atoms '+str((a1,a2))+' ... bond_length, batom_names = '+str(return_val)+'\n')
    return return_val








def LookupRestAngle(a1, a2, a3,
                    atom2equiv_angle,
                    angle2theta0,
                    atom2auto_angle,
                    angle2theta0_auto):
    """ 
    Try to find angle parameters between atoms whose original atom
    names (without equivalences) are a1, a2, and a3.  Then return
    both the equilibrium rest angle for that 3body interaction
    as well as the equivalent atom names used to look it up. (These
    could be stored in either atom2equiv_angle or atom2auto_angle.)
    If a match was not found, return None.
    """
    return_val = None
    anames = (atom2equiv_angle[a1], atom2equiv_angle[a2], atom2equiv_angle[a3])
    angle_name = EncodeInteractionName(ReverseIfEnds(anames))
    if angle_name in angle2theta0:
        return_val = (angle2theta0[angle_name], [anames[0], anames[1], anames[2]])

    # If no angle between these atoms is defined, 
    # check the angles in the _auto section(s)
    # This is a lot messier.
    elif ((a1 in atom2auto_angle[0]) and
          (a2 in atom2auto_angle[1]) and
          (a3 in atom2auto_angle[2])):

        anames = [atom2auto_angle[0][a1],
                  atom2auto_angle[1][a2],
                  atom2auto_angle[2][a3]]
        #sys.stderr.write('DEBUG: LookupRestAngle(): a1,a2,a3=('+
        #                 a1+','+a2+','+a3+'), anames='+str(anames)+'\n')

        # Because _auto interactions can contain wildcards,
        # there can be multiple entries in angle2theta0_auto[]
        # for the same list of atom names, and we have to
        # consider all of them, and pick the one with the
        # most priority (ie. whose priority number is lowest).
        # (Note: The MSI file format uses low priority numbers
        #  to indicate high priority.  Somewhat confusing.)
        HUGE_VAL = 2000000000
        best_priority = HUGE_VAL
        pattern = ['','','']
        for (pattern[0],pattern[1],pattern[2]), theta0 in angle2theta0_auto.items():
            priority = DoAtomsMatchPattern(anames, pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (theta0, [anames[0], anames[1], anames[2]])
            anames.reverse() # now check of the atoms in reverse order match
            priority = DoAtomsMatchPattern(anames, pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (theta0, [anames[2], anames[1], anames[0]]) #preserve atom order
        #if return_val != None:
        #    sys.stderr.write('DEBUG: For atoms '+str((a1,a2))+' ... rest_angle, anames = '+str(return_val)+'\n')
    return return_val


                              





def Equivalences2ffids(lines_equivalences,
                       atom_types,
                       atom2equiv_pair,
                       atom2equiv_bond,
                       atom2equiv_angle,
                       atom2equiv_dihedral,
                       atom2equiv_improper):
    """
    This function reads a list of lines containing "equivalences" and
    "auto_equivalences" from an MSI-formatted .FRC file.
    Then, for each atom type, it generates a long string which includes the 
    original atom type name as well as all of the equivalences it belongs to.
    Later on, when it is time to generate angles, dihedrals, or impropers,
    moltemplate will search for patterns contained in these strings to decide
    which type of interaction to generate.
    This function returns a dictionary that converts the original atom type name
    into these strings.
    """
    for line in lines_equivalences:
        #tokens = SplitQuotedString(line.strip(),
        #                           comment_char='!>')

        # skip past both '!' and '>' characters
        ic1 = line.find('!')
        ic = ic1
        ic2 = line.find('>')
        if ic2 != -1 and ic2 < ic1:
            ic = ic2
        if ic != -1:
            line = line[:ic]
        else:
            line = line.rstrip('\n')
        tokens = line.strip().split()
        #sys.stderr.write('DEBUG Equivalences2ffids():\n'
        #                 '      tokens = '+str(tokens)+'\n')
        atype = tokens[2]
        atom2equiv_pair[atype] = tokens[3]
        atom2equiv_bond[atype] = tokens[4]
        atom2equiv_angle[atype] = tokens[5]
        atom2equiv_dihedral[atype] = tokens[6]
        atom2equiv_improper[atype] = tokens[7]

    atom2ffid = OrderedDict()
    for atom in atom_types:
        atom2ffid[atom] = (atom + 
                           #',p'+atom2equiv_pair.get(atom,'') + 
                           ',b'+atom2equiv_bond.get(atom,'') + 
                           ',a'+atom2equiv_angle.get(atom,'') + 
                           ',d'+atom2equiv_dihedral.get(atom,'') + 
                           ',i'+atom2equiv_improper.get(atom,''))
    return atom2ffid






def AutoEquivalences2ffids(lines_equivalences,
                           lines_auto_equivalences,
                           atom_types,
                           atom2equiv_pair,
                           atom2equiv_bond,
                           atom2equiv_angle,
                           atom2equiv_dihedral,
                           atom2equiv_improper,
                           atom2auto_pair,
                           atom2auto_bondincr,
                           atom2auto_bond,
                           atom2auto_angleend,
                           atom2auto_anglecenter,
                           atom2auto_dihedralend,
                           atom2auto_dihedralcenter,
                           atom2auto_improperend,
                           atom2auto_impropercenter):
    """
    This function is a variant of Equivalences2ffids() which also considers
    "auto_equivalences".
    This function returns a dictionary that converts the original atom type name
    into a string that includes that atom's "equivalences",
    as well as its "auto_equivalences".
    moltemplate will search for patterns contained in these strings to decide
    which type of interaction to generate.
    """
    Equivalences2ffids(lines_equivalences,
                       atom_types,
                       atom2equiv_pair,
                       atom2equiv_bond,
                       atom2equiv_angle,
                       atom2equiv_dihedral,
                       atom2equiv_improper)

    # ------ The following lines are for processing "auto_equivalences" -----
    #
    # What is the difference between "equivalences" and "auto_equivalences"?
    #
    # equivalences:
    # Here is an excerpt from the Discover manual describing "equivalences":
    #  "Chemically distinct atoms often differ in some, but not all,
    #   of their forcefield parameters. For example, the bond parameters
    #  for the C-C bonds in ethene and in benzene are quite different,
    #  but the nonbond parameters for the carbon atoms are essentially
    #  the same. Rather than duplicating the nonbond parameters in the
    #  forcefield parameter file, the Discover program uses atom type
    #  equivalences to simplify the problem. In the example, the phenyl
    #  carbon atom type is equivalent to the pure sp2 carbons of ethene
    #  insofar as the nonbond parameters are concerned. The Discover
    #  program recognizes five types of equivalences for each atom
    #  type: nonbond, bond, angle, torsion, and out-of-plane.
    #  Cross terms such as bond-bond terms have the same equivalences
    #  (insofar as atom types are concerned) as the diagonal term of
    #  the topology of all the atoms defining the internal coordinates.
    #  For the bond-bond term, this means that the atom type
    #  equivalences for angles would be used
    #
    # auto_equivalences:
    #   Are similar to equivalences, but apparently with lower priority.
    #   In addition, it seems that, when looking up some of the class2 terms
    #   in the interaction according to atom type using "auto_equivalences"
    #   a distinction is made between end atoms and central atoms.
    #   The parameters for these interactions are also stored in different 
    #   tables in the .frc file, with different comments/tags.
    #   (for example, "cff91_auto" as opposed to "cff91")
    # An excerpt from the Discover manual is somewhat vague:
    #  "A forcefield may include automatic parameters for use when
    #   better-quality explicit parameters are not defined for a
    #   particular bond, angle, torsion, or out-of-plane interaction.
    #   These parameters are intended as temporary patches, to allow
    #   you to begin calculations immediately."

    for line in lines_auto_equivalences:
        #tokens = SplitQuotedString(line.strip(),
        #                           comment_char='!>')

        # skip past both '!' and '>' characters
        ic1 = line.find('!')
        ic = ic1
        ic2 = line.find('>')
        if ic2 != -1 and ic2 < ic1:
            ic = ic2
        if ic != -1:
            line = line[:ic]
        else:
            line = line.rstrip('\n')
        tokens = line.strip().split()
        #sys.stderr.write('DEBUG Equivalences2ffids():\n'
        #                 '      tokens = '+str(tokens)+'\n')
        atype = tokens[2]
        atom2auto_pair[atype] = tokens[3]
        atom2auto_bondincr[atype] = tokens[4]
        atom2auto_bond[atype] = tokens[5]
        atom2auto_angleend[atype] = tokens[6]
        atom2auto_anglecenter[atype] = tokens[7]
        atom2auto_dihedralend[atype] = tokens[8]
        atom2auto_dihedralcenter[atype] = tokens[9]
        atom2auto_improperend[atype] = tokens[10]
        atom2auto_impropercenter[atype] = tokens[11]

    atom2ffid = OrderedDict()
    for atom in atom_types:
        atom2ffid[atom] = (atom + 
                           #',p'+atom2equiv_pair.get(atom,'') + 
                           ',b'+atom2equiv_bond.get(atom,'') + 
                           ',a'+atom2equiv_angle.get(atom,'') + 
                           ',d'+atom2equiv_dihedral.get(atom,'') + 
                           ',i'+atom2equiv_improper.get(atom,'') + 
                           #',ap'+atom2auto_pair.get(atom,'') + 
                           #',aq'+atom2auto_bondincr.get(atom,'') + 
                           #',ab'+atom2auto_bond.get(atom,'') + 
                           ',aae'+atom2auto_angleend.get(atom,'') + 
                           ',aac'+atom2auto_anglecenter.get(atom,'') + 
                           ',ade'+atom2auto_dihedralend.get(atom,'') + 
                           ',adc'+atom2auto_dihedralcenter.get(atom,'') + 
                           ',aie'+atom2auto_improperend.get(atom,'') + 
                           ',aic'+atom2auto_impropercenter.get(atom,'') +
                           ''
                          )
    return atom2ffid






def main():
    try:
        sys.stderr.write(g_program_name + ", version " +
                         __version__ + ", " + __date__ + "\n")
        if sys.version < '2.6':
            raise InputError('Error: Using python ' + sys.version + '\n' +
                             '       Alas, your version of python is too old.\n'
                             '       You must upgrade to a newer version of python (2.6 or later).')
    
        if sys.version < '2.7':
            from ordereddict import OrderedDict
        else:
            from collections import OrderedDict 
    
        if sys.version > '3':
            import io
        else:
            import cStringIO
    
        # defaults:
        ffname = 'BIOSYM_MSI_FORCE_FIELD'
        type_subset = set([])
        filename_in = ''
        file_in = sys.stdin
        #file_in = open('compass_published.frc','r')  #CONTINUEHERE
        include_auto_equivalences = True
        #pair_style_name = 'lj/class2/coul/long'
        #pair_style_params = "10.0 10.0"
        pair_style2docs = {}
        pair_style2args = defaultdict(str)
        pair_style2docs['lj/cut/coul/long'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pair_style2args['lj/cut/coul/long'] = '10.0'
        pair_style2docs['lj/class2/coul/long'] = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pair_style2args['lj/class2/coul/long'] = '10.0'
        pair_style2docs['lj/class2/coul/cut'] = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pair_style2args['lj/class2/coul/cut'] = '10.0'

        bond_style2docs = {}
        #bond_style2args = defaultdict(str)
        bond_style2docs['harmonic'] = 'http://lammps.sandia.gov/doc/bond_harmonic.html'
        bond_style2docs['class2'] = 'http://lammps.sandia.gov/doc/bond_class2.html'
        bond_style2docs['morse'] = 'http://lammps.sandia.gov/doc/bond_morse.html'

        angle_style2docs = {}
        #angle_style2args = defaultdict(str)
        angle_style2docs['harmonic'] = 'http://lammps.sandia.gov/doc/angle_harmonic.html'
        angle_style2docs['class2'] = 'http://lammps.sandia.gov/doc/angle_class2.html'

        dihedral_style2docs = {}
        #dihedral_style2args = defaultdict(str)
        dihedral_style2docs['charmm'] = 'http://lammps.sandia.gov/doc/dihedral_charmm.html'
        dihedral_style2docs['class2'] = 'http://lammps.sandia.gov/doc/dihedral_class2.html'
        dihedral_symmetry_subgraph = ''  # default

        improper_style2docs = {}
        #improper_style2args = defaultdict(str)
        improper_style2docs['cvff'] = 'http://lammps.sandia.gov/doc/improper_cvff.html'
        improper_style2docs['class2'] = 'http://lammps.sandia.gov/doc/improper_class2.html'
        improper_symmetry_subgraph = 'cenJsortIKL'

        pair_mixing_style = 'sixthpower tail yes'

        special_bonds_command = 'special_bonds lj/coul 0.0 0.0 1.0 dihedral yes'
        # Thanks to Paul Saxe for is suggestions
        # http://lammps.sandia.gov/threads/msg11270.html


        kspace_style = 'kspace_style pppm 0.0001'
        pair_styles_selected = set([])
        #pair_style_link = 'http://lammps.sandia.gov/doc/pair_class2.html'
        #pair_style_args = '10.0'
        #pair_style_command = "    pair_style hybrid " + \
        #    pair_style_name + " " + pair_style_args + "\n"
        bond_styles_selected = set([])
        #bond_style_link = bond_style2docs[bond_style_name]
        #bond_style_args = ''
        angle_styles_selected = set([])
        #angle_style_link = angle_style2docs[angle_style_name]
        #angle_style_args = ''
        dihedral_styles_selected = set([])
        #dihedral_style_link = dihedral_style2docs[dihedral_style_name]
        #dihedral_style_args = ''
        improper_styles_selected = set([])
        #improper_style_link = improper_style2docs[improper_style_name]
        #improper_style_args = ''
        hbond_style_name = ''
        hbond_style_link = ''
        hbond_style_args = ''

        
    
        argv = [arg for arg in sys.argv]
    
        i = 1
    
        while i < len(argv):
    
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
    
            if argv[i] == '-atoms':
                if i + 1 >= len(argv):
                    raise InputError('Error: the \"' + argv[i] + '\" argument should be followed by a quoted string\n'
                                     '       which contains a space-delimited list of of a subset of atom types\n'
                                     '       you want to use from the original force-field.\n'
                                     '       Make sure you enclose the entire list in quotes.\n')
                type_subset = set(argv[i + 1].strip('\"\'').strip().split())
                del argv[i:i + 2]
    
            elif argv[i] == '-name':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of the force-field\n')
                ffname = argv[i + 1]
                del argv[i:i + 2]
    
            elif argv[i] in ('-file', '-in-file'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of a force-field file\n')
                filename_in = argv[i + 1]
                try:
                    file_in = open(filename_in, 'r')
                except IOError:
                    sys.stderr.write('Error: Unable to open file\n'
                                     '       \"' + filename_in + '\"\n'
                                     '       for reading.\n')
                    sys.exit(1)
                del argv[i:i + 2]
    
            elif argv[i] == '-pair-cutoff':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a number'
                                     '       (or two numbers enclosed in a single pair of quotes)\n')
                pair_style_args = argv[i+1]
                del argv[i:i + 2]

            elif argv[i] == '-pair-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by either \"lj/class2/coul/cut\" or \"lj/class2/coul/long\"\n')
                pair_styles = argv[i + 1].split(',')
                for pair_style in pair_styles:
                    if pair_style == '9-6':
                        pair_style = 'lj/class2/coul/long'
                    elif pair_style in ('12-6', 'lj', 'LJ'):
                        pair_style = 'lj/cut/coul/long'

                    if  pair_style.find('lj/class2/coul/long') == 0:
                        kspace_style = 'kspace_style pppm 0.0001'
                    elif pair_style.find('lj/cut/coul/long') == 0:
                        kspace_style = 'kspace_style pppm 0.0001'
                    elif pair_style.find('lj/class2/coul/cut') == 0:
                        pass
                        #kspace_style = ''
                    elif pair_style.find('lj/cut') == 0:
                        pass
                        #kspace_style = ''
                    else:
                        raise InputError('Error: ' + argv[i] + ' ' + pair_style_name + ' not supported.\n'
                                         '          The following pair_styles are supported:\n'
                                         '       lj/class2/coul/cut\n'
                                         '       lj/class2/coul/long\n'
                                         '       lj/cut\n'
                                         '       lj/cut/coul/long\n')
                    pair_styles_selected.add(pair_style)

                del argv[i:i + 2]

            elif argv[i] == '-bond-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible bond_style.\n')
                bond_styles = argv[i + 1].split(',')
                for bond_style in bond_styles:
                    bond_styles_selected.add(bond_style)
                #bond_style2args[bond_style] = argv[i + 1].split()[1:]
                #if bond_style_name.find('harmonic') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_harmonic.html'
                #elif bond_style_name.find('morse') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_morse.html'
                #elif bond_style_name.find('class2') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\", \"class2\", or \"morse\".\n')
                del argv[i:i + 2]

            elif argv[i] == '-angle-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible angle_style.\n')
                angle_styles = argv[i + 1].split(',')
                for angle_style in angle_styles:
                    angle_styles_selected.add(angle_style)
                #if angle_style_name.find('harmonic') == 0:
                #    pass
                #    #angle_style_link = 'http://lammps.sandia.gov/doc/angle_harmonic.html'
                #elif angle_style_name.find('class2') == 0:
                #    pass
                #    #angle_style_link = 'http://lammps.sandia.gov/doc/angle_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-dihedral-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible dihedral_style.\n')
                dihedral_styles = argv[i + 1].split(',')
                for dihedral_style in dihedral_styles:
                    dihedral_styles_selected.add(dihedral_style)
                #if dihedral_style_name.find('charmm') == 0:
                #    pass
                #    #dihedral_style_link = 'http://lammps.sandia.gov/doc/dihedral_charmm.html'
                #elif dihedral_style_name.find('class2') == 0:
                #    pass
                #    #dihedral_style_link = 'http://lammps.sandia.gov/doc/dihedral_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-improper-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible impropoer_style.\n')
                improper_styles = argv[i + 1].split(',')
                for improper_style in improper_styles:
                    improper_styles_selected.add(improper_style)
                #if impropoer_style_name.find('harmonic') == 0:
                #    pass
                #    #impropoer_style_link = 'http://lammps.sandia.gov/doc/impropoer_harmonic.html'
                #elif impropoer_style_name.find('class2') == 0:
                #    pass
                #    #impropoer_style_link = 'http://lammps.sandia.gov/doc/impropoer_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-hbond-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' ' + hbond_style_name + '\n'
                                     '       should be followed by a compatible pair_style.\n')
                hbond_style_name = argv[i + 1]
                hbond_style_link = 'http://lammps.sandia.gov/doc/pair_hbond_dreiding.html'
                if hbond_style_name.find('none') == 0:
                    hbond_style_name = ''
                    hbond_style_args = ''
                elif hbond_style_name.find('hbond/dreiding/lj') == 0:
                    n = len('hbond/dreiding/lj')
                    hbond_style_args = hbond_style_name[n+1:]
                    hbond_style_name = hbond_style_name[:n]
                elif hbond_style_name.find('hbond/dreiding/morse') == 0:
                    n = len('hbond/dreiding/morse')
                    hbond_style_args = hbond_style_name[n+1:]
                    hbond_style_name = hbond_style_name[:n]
                else:
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by either\n'
                                     '       \"hbond/dreiding/lj\" or \"hbond/dreiding/morse"\n')
                del argv[i:i + 2]

            elif argv[i] in ('-url', '-in-url'):
                import urllib2
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of a\n'
                                     '       file containing force-field information in msi/frc format.\n')
                url = argv[i + 1]
                try:
                    request = urllib2.Request(url)
                    file_in = urllib2.urlopen(request)
                except urllib2.URLError:
                    sys.stdout.write("Error: Unable to open link:\n" + url + "\n")
                    sys.exit(1)
                del argv[i:i + 2]
    
            #elif argv[i] == '-auto':
            #    include_auto_equivalences = True
            #    del argv[i:i + 1]
    
            elif argv[i] in ('-help', '--help', '-?', '--?'):
                sys.stderr.write(doc_msg)
                sys.exit(0)
                del argv[i:i + 1]
    
            else:
                i += 1
    
        if len(argv) != 1:
            raise InputError('Error: Unrecongized arguments: ' + ' '.join(argv[1:]) +
                             '\n\n' + doc_msg)

        # Default styles:
        if len(bond_styles_selected) == 0:
            bond_styles_selected.add('class2')
        if len(angle_styles_selected) == 0:
            angle_styles_selected.add('class2')
        if len(dihedral_styles_selected) == 0:
            dihedral_styles_selected.add('class2')
        if len(improper_styles_selected) == 0:
            improper_styles_selected.add('class2')
        if len(pair_styles_selected) == 0:
            pair_styles_selected.add('lj/class2/coul/long')

        #sys.stderr.write("Reading parameter file...\n")

        lines = file_in.readlines()
        atom2charge = OrderedDict()  # lookup charge from atom type
        atom2mass = OrderedDict()  # lookup mass from atom type
        # equivalences lookup
        atom2ffid = OrderedDict()  # lookup "force-field-ID" a string containing
                                   # equivalences to lookup bonded interactions
        atom2equiv_pair = OrderedDict() # lookup the equivalent symbol used for
                                        # looking up pair interactions
        atom2equiv_bond = OrderedDict()
        atom2equiv_angle = OrderedDict()
        atom2equiv_dihedral = OrderedDict()
        atom2equiv_improper = OrderedDict()
        # inverse equivalences lookup
        equiv_pair2atom = defaultdict(set)
        equiv_bond2atom = defaultdict(set)
        equiv_angle2atom = defaultdict(set)
        equiv_dihedral2atom = defaultdict(set)
        equiv_improper2atom = defaultdict(set)

        # auto equivalences lookup
        atom2auto_pair = OrderedDict()
        atom2auto_bondincr = OrderedDict()
        atom2auto_bond = OrderedDict()
        atom2auto_angleend = OrderedDict()
        atom2auto_anglecenter = OrderedDict()
        atom2auto_dihedralend = OrderedDict()
        atom2auto_dihedralcenter = OrderedDict()
        atom2auto_improperend = OrderedDict()
        atom2auto_impropercenter = OrderedDict()
        # inverse auto equivalences lookup
        auto_pair2atom = defaultdict(set)
        auto_bondincr2atom = defaultdict(set)
        auto_bond2atom = defaultdict(set)
        auto_angleend2atom = defaultdict(set)
        auto_anglecenter2atom = defaultdict(set)
        auto_dihedralend2atom = defaultdict(set)
        auto_dihedralcenter2atom = defaultdict(set)
        auto_improperend2atom = defaultdict(set)
        auto_impropercenter2atom = defaultdict(set)


        atom2element = OrderedDict()  # Optional:
                                      # which element (eg 'C', 'O') ? (Note this
                                      # is different from atom type: 'C1', 'Oh')
        atom2numbonds = OrderedDict() # Optional: how many bonds emanate from
        atom2descr = OrderedDict()    # Optional: a brief description
        atom2ver = OrderedDict()  # atoms introduced in different versions of ff
        atom2ref = OrderedDict()  # reference to paper where atom introduced
        lines_equivalences = []      # equivalences for force-field lookup
        lines_auto_equivalences = [] # auto_equivalences have lower priority

        pair2params = OrderedDict()
        pair2style = OrderedDict()
        pair_styles = set([])
        pair2ver = OrderedDict()
        pair2ref = OrderedDict()

        bond2chargepair = OrderedDict()      # a.k.a "bond increments"
        charge_pair_priority = OrderedDict() # priority in case multiple entries
                                             # exist for the same pair of atoms
        charge_pair_ver = OrderedDict        # which version of the force field?
        charge_pair_ref = OrderedDict        # paper introducing this chargepair

        bond2params = OrderedDict()  # store a tuple with the 2-body bond
                                     # interaction type, and its parameters
                                     # for every type of bond
        bond2priority = OrderedDict() # What is the priority of this interaction?
        bond2style = OrderedDict()    # What LAMMPS bond style (formula)
                                      # is used for a given interaction?
        bond_styles = set([])         # Contains all bond styles used.
        bond2ver = OrderedDict()
        bond2ref = OrderedDict()
        bond2r0 = OrderedDict()
        bond2r0_auto = OrderedDict()

        angle2params = OrderedDict() # store a tuple with the 3-body angle
                                     # interaction type, and its parameters
                                     # for every type of angle

        # http://lammps.sandia.gov/doc/angle_class2.html
        #angle2class2_a = OrderedDict()  # params for the "a" class2 terms
        angle2class2_bb = OrderedDict() # params for the "bb" class2 terms
        angle2class2_bb_sh = OrderedDict()
        angle2class2_ba = OrderedDict() # params for the "ba" class2 terms
        angle2class2_ba_sh = OrderedDict()
        angle2priority = OrderedDict()  # What is the priority of this interaction?
        angle2priority_sh = OrderedDict()
        angle2style = OrderedDict()    # What LAMMPS angle style (formula)
                                       # is used for a given interaction?
        angle_styles = set([])         # Contains all angle styles used.
        angle2ver = OrderedDict()
        angle2ver_sh = OrderedDict()
        angle2ver_bb_sh = OrderedDict()
        angle2ver_ba_sh = OrderedDict()
        angle2ref = OrderedDict()
        angle2theta0 = OrderedDict()
        angle2theta0_auto = OrderedDict()

        # http://lammps.sandia.gov/doc/dihedral_class2.html
        dihedral2params = OrderedDict() # store a tuple with the 4-body dihedral
                                        # interaction type, and its parameters
                                        # for every type of dihedral
        dihedral2params_sh = OrderedDict()
        #dihedral2class2_d = OrderedDict() # params for the "d" class2 term
        dihedral2class2_mbt = OrderedDict() # params for the "mbt" class2 term
        dihedral2class2_mbt_sh = OrderedDict()
        dihedral2class2_ebt = OrderedDict() # params for the "ebt" class2 term
        dihedral2class2_ebt_sh = OrderedDict()
        #dihedral2sym_ebt = OrderedDict()
        dihedral2class2_at = OrderedDict() # params for the "at" class2 term
        dihedral2class2_at_sh = OrderedDict()
        #dihedral2sym_at = OrderedDict()
        dihedral2class2_aat = OrderedDict() # params for the "aat" class2 term
        dihedral2class2_aat_sh = OrderedDict()
        #dihedral2sym_aat = OrderedDict()
        dihedral2class2_bb13 = OrderedDict() # params for the "bb13" class2 term
        dihedral2class2_bb13_sh = OrderedDict()
        #dihedral2sym_bb13 = OrderedDict()
        dihedral2priority = OrderedDict()  # What is the priority of this interaction?
        dihedral2priority_sh = OrderedDict()
        dihedral2style = OrderedDict()    # What LAMMPS dihedral style (formula)
                                          # is used for a given interaction?
        dihedral_styles = set([])         # Contains all dihedral styles used.
        dihedral2ver_sh = OrderedDict()
        dihedral2ver_mbt_sh = OrderedDict()
        dihedral2ver_ebt_sh = OrderedDict()
        dihedral2ver_at_sh = OrderedDict()
        dihedral2ver_aat_sh = OrderedDict()
        dihedral2ver_bb13_sh = OrderedDict()
        dihedral2ref = OrderedDict()


        # http://lammps.sandia.gov/doc/improper_class2.html
        improper2params = OrderedDict() # store a tuple with the 4-body improper
                                        # interaction type, and its parameters
                                        # for every type of imporpoer
        #improper2class2_i = OrderedDict() # params for the "i" class2 term
        improper2class2_aa = OrderedDict() # params for the "aa" class2 term

        improper2cross = defaultdict(dict)
                           # improper2cross[imp_name][atoms] stores the 
                           # coefficient (K) for the angle-angle ("aa") 
                           # improper interactions between a pair of 
                           # neighboring 3-body angles (in the .FRC file).
                           # "imp_name" is the name of the improper interaction
                           #   (which is a concatination of the central atom and
                           #   the 3 surrounding leaf atoms (which are sorted))
                           # "atoms" indicates, for that K value, the list of
                           #   leaf atoms for that K value as they appear in the
                           #   corresponding line of the .frc file (however the
                           #   and last atom names are swapped if the first
                           #   atom name is lexicographically > the last, to
                           #   eliminate redundancy and ambiguity.)

        improper2sym = defaultdict(set)
                           # improper2sym[imp_name] indicates which subset of
                           # leaf atoms (from 0 to 2) are equivalent and can
                           # tolerate having their order rearranged without
                           # effecting the energy.  Later on this will be used
                           # to reduce the number of improper interactions that
                           # will be generated by moltemplate.

        improper2priority = OrderedDict() # What is the priority of this interaction?
        improper2priority_sh = OrderedDict()
        improper2style = OrderedDict()    # What LAMMPS improper style (formula)
                                          # is used for a given interaction?
        improper_styles = set([])         # Contains all improper styles used.
        improper2ver_sh = OrderedDict()
        improper2ver_aa = OrderedDict()        
        improper2ver = OrderedDict()
        improper2ref = OrderedDict()


        # Warn users if force field contains terms which cannot yet
        # be simulated with LAMMPS (as of 2017-2-07)
        display_OOP_OOP_warning = False
        display_torsion_torsion_1_warning = False


        """
         --- these next few lines of code appear to be unnecessary.
         --- I'll probably delete this code in a later version
        hbond2params = OrderedDict()    # lookup hbond parameters and atom types
        hbond2donors = OrderedDict()    # according to the identifier in the 2nd
        hbond2acceptors = OrderedDict() #  column of the "#hbond_definition"
        hbond2hydrogens = OrderedDict() # section of an .frc file.
        """

        allowed_section_names = set(['#define',
                                     # sections used in all MSI force-fields
                                     '#atom_types',
                                     '#equivalence',
                                     '#auto_equivalence',
                                     '#nonbond(9-6)',
                                     '#nonbond(12-6)',
                                     '#quadratic_bond',
                                     '#quartic_bond',
                                     '#morse_bond',
                                     '#quadratic_angle',
                                     '#quartic_angle',
                                     '#bond-bond',
                                     '#bond-angle',
                                     '#torsion_1',
                                     '#torsion_3',
                                     '#middle_bond-torsion_3',
                                     '#end_bond-torsion_3',
                                     '#angle-torsion_3',
                                     '#angle-angle-torsion_1',#(class2 dihedral)
                                     '#bond-bond_1_3', #(a class2 dihedral term)
                                     '#out_of_plane',
                                     '#wilson_out_of_plane',
                                     '#angle-angle',   #(a class2 improper term)
                                     '#out_of_plane-out_of_plane', # UNSUPPORTED
                                     '#torsion-torsion_1',         # UNSUPPORTED
                                     '#bond_increments',
                                     '#hbond_definition',          # irrelevant?
                                     '#templates',
                                     '#reference',
                                     '#end'
                                     ])

        icol_type = icol_mass = icol_elem = icol_nbonds = icol_comment = icol_ver = icol_ref = -1

        section_name = ''
        section_is_auto = False

        sys.stderr.write("parsing file pass1: look for atom types and equivalences...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if tokens[0] in allowed_section_names:
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif not tokens[0] in ('#version',
                                       '#define'):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')
            elif (len(tokens) == 8) and (section_name == '#equivalence'):
                if line.lstrip().find('!') == 0:
                    continue
                lines_equivalences.append(line)
            elif (len(tokens) == 12) and (section_name == '#auto_equivalence'):
                if line.lstrip().find('!') == 0:
                    continue
                lines_auto_equivalences.append(line)
            elif (len(tokens) > 0) and (section_name == '#atom_types'):
                # Different FRC files put this information in different
                # columns.  Column order is stored in the !Ver comment line:
                if line.lstrip().find('!Ver') == 0:
                    tokens = line.strip().split()
                    for i in range(0, len(tokens)):
                        if tokens[i].lower() == 'type':
                            icol_type = i
                        elif tokens[i].lower() == 'mass':
                            icol_mass = i
                        elif tokens[i].lower() == 'element':
                            icol_elem = i
                        elif tokens[i].lower() == 'connections':
                            icol_nbonds = i
                        elif tokens[i].lower() == 'comment':
                            icol_comment = i
                        elif tokens[i].lower() == '!ver':   #(version of ff)
                            icol_ver = i
                        elif tokens[i].lower() == 'ref':
                            icol_ref = i
                    assert(icol_ver == 0)

                    if -1 in (icol_type, icol_mass):
                        raise InputError('Error: Invalid #atom_types section.\n'
                                         '       The meaning of each column cannot be determined.\n'
                                         '       This file needs a valid "!Ver..." comment.\n')
                    if icol_comment == -1:
                        icol_comment = max(icol_type, icol_mass,
                                           icol_elem, icol_nbonds) + 1

                    sys.stderr.write('icol_ver = '+str(icol_ver)+'\n')
                    sys.stderr.write('icol_ref = '+str(icol_ref)+'\n')
                    sys.stderr.write('icol_mass = '+str(icol_mass)+'\n')
                    sys.stderr.write('icol_nelem = '+str(icol_elem)+'\n')
                    sys.stderr.write('icol_nbonds = '+str(icol_nbonds)+'\n')
                    sys.stderr.write('icol_comment = '+str(icol_comment)+'\n')
                    continue

                tokens = map(RemoveOuterQuotes,
                             NSplitQuotedString(line.strip(),
                                                icol_comment+1,
                                                quotes='',
                                                comment_char='>'))
                tokens = list(tokens)

                if (len(tokens) > 4):
                    if ((len(type_subset) == 0) or (tokens[1] in type_subset)):

                        atom2mass[tokens[icol_type]] = str(max(float(tokens[icol_mass]), 1.0e-06))
                        # Some atoms in cvff.prm have zero mass. Unfortunately this
                        # causes LAMMPS to crash, even if these atoms are never used,
                        # so I give the mass a non-zero value instead.

                        if icol_elem != -1:
                            atom2element[tokens[icol_type]] = tokens[icol_elem]
                        if icol_nbonds != -1:
                            atom2numbonds[tokens[icol_type]] = int(tokens[icol_nbonds])
                        atom2descr[tokens[icol_type]] = tokens[icol_comment]
                        atom2ver[tokens[icol_type]] = tokens[icol_ver]
                        atom2ref[tokens[icol_type]] = tokens[icol_ref]

                elif len(tokens) > 0:
                    raise InputError('Error: Invalid atom line: (line#'+str(iline)+')\n' +
                                     '\"'+line.strip()+'\"')

        atom_types = [x for x in atom2mass]

        # Now construct the lookup tables and inverse tables
        # we will need to understand the remainder of the file:
        if not include_auto_equivalences:
            atom2ffid = Equivalences2ffids(lines_equivalences,
                                           atom_types,
                                           atom2equiv_pair,
                                           atom2equiv_bond,
                                           atom2equiv_angle,
                                           atom2equiv_dihedral,
                                           atom2equiv_improper)
        else:
            atom2ffid = AutoEquivalences2ffids(lines_equivalences,
                                               lines_auto_equivalences,
                                               atom_types,
                                               atom2equiv_pair,
                                               atom2equiv_bond,
                                               atom2equiv_angle,
                                               atom2equiv_dihedral,
                                               atom2equiv_improper,
                                               atom2auto_pair,
                                               atom2auto_bondincr,
                                               atom2auto_bond,
                                               atom2auto_angleend,
                                               atom2auto_anglecenter,
                                               atom2auto_dihedralend,
                                               atom2auto_dihedralcenter,
                                               atom2auto_improperend,
                                               atom2auto_impropercenter)

        for a,e in atom2equiv_pair.items():
            equiv_pair2atom[e].add(a)
        for a,e in atom2equiv_bond.items():
            equiv_bond2atom[e].add(a)
        for a,e in atom2equiv_angle.items():
            equiv_angle2atom[e].add(a)
        for a,e in atom2equiv_dihedral.items():
            equiv_dihedral2atom[e].add(a)
        for a,e in atom2equiv_improper.items():
            equiv_improper2atom[e].add(a)

        # the inverse lookup for '*' matches all atom types
        for a in atom_types:
            equiv_pair2atom['*'].add(a)
            equiv_pair2atom['X'].add(a)
            equiv_bond2atom['*'].add(a)
            equiv_bond2atom['X'].add(a)
            equiv_angle2atom['*'].add(a)
            equiv_angle2atom['X'].add(a)
            equiv_dihedral2atom['*'].add(a)
            equiv_dihedral2atom['X'].add(a)
            equiv_improper2atom['*'].add(a)
            equiv_improper2atom['X'].add(a)

        for a,e in atom2auto_pair.items():
            auto_pair2atom[e].add(a)
        for a,e in atom2auto_bondincr.items():
            auto_bondincr2atom[e].add(a)
        for a,e in atom2auto_bond.items():
            auto_bond2atom[e].add(a)
        for a,e in atom2auto_angleend.items():
            auto_angleend2atom[e].add(a)
            #auto_angle[0][e].add(a)
            #auto_angle[2][e].add(a)
        for a,e in atom2auto_anglecenter.items():
            auto_anglecenter2atom[e].add(a)
            #auto_angle[1][e].add(a)
        for a,e in atom2auto_dihedralend.items():
            auto_dihedralend2atom[e].add(a)
            #auto_dihedral2atom[0][e].add(a)
            #auto_dihedral2atom[3][e].add(a)
        for a,e in atom2auto_dihedralcenter.items():
            auto_dihedralcenter2atom[e].add(a)
            #auto_dihedral2atom[1][e].add(a)
            #auto_dihedral2atom[2][e].add(a)
        for a,e in atom2auto_improperend.items():
            auto_improperend2atom[e].add(a)
        for a,e in atom2auto_impropercenter.items():
            auto_impropercenter2atom[e].add(a)

        # the inverse lookup for '*' matches all atom types
        for a in atom_types:
            auto_pair2atom['*'].add(a)
            auto_pair2atom['X'].add(a)
            auto_bondincr2atom['*'].add(a)
            auto_bondincr2atom['X'].add(a)
            auto_bond2atom['*'].add(a)
            auto_bond2atom['X'].add(a)
            auto_angleend2atom['*'].add(a)
            auto_angleend2atom['X'].add(a)
            auto_anglecenter2atom['*'].add(a)
            auto_anglecenter2atom['X'].add(a)
            auto_dihedralend2atom['*'].add(a)
            auto_dihedralend2atom['X'].add(a)
            auto_dihedralcenter2atom['*'].add(a)
            auto_dihedralcenter2atom['X'].add(a)
            auto_improperend2atom['*'].add(a)
            auto_improperend2atom['X'].add(a)
            auto_impropercenter2atom['*'].add(a)
            auto_impropercenter2atom['X'].add(a)








        sys.stderr.write("parsing file pass2: look for bonds, bond_increments and nonbonded (pair) interactions...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')


            elif ((len(tokens) > 4) and (section_name == '#nonbond(12-6)')
                  and (pair_styles_selected & set(['lj','lj/cut','lj/cut/coul/long',
                                                   'lj/cut/coul/cut','lj/cut/coul/debye',
                                                   'lj/cut/coul/dsf','lj/cut/coul/msm',
                                                   '12-6','nonbond(12-6)']))):

                if line.lstrip().find('!') == 0:
                    continue
                atom_name = tokens[2]
                pair2ver[atom_name] = tokens[0]
                pair2ref[atom_name] = tokens[1]
                A = float(tokens[3])
                B = float(tokens[4])
                epsilon = B*B/(4*A)
                sigma = pow(B/A, 1.0/6)
                if sigma == 0.0:
                    sigma = 1.0   #(non-zero to avoid nan error later)
                pair2style[atom_name] = 'lj/cut/coul/long'
                pair_styles.add('lj/cut/coul/long')
                pair2params[atom_name] = (str(epsilon)+' '+str(sigma))
                pair_mixing_style = 'geometric tail yes'
                #if pair_style_name.find('lj/cut') == 0:
                #    pair2params[atom_name] = (str(epsilon)+' '+str(sigma))
                #    pair_mixing_style = 'geometric tail yes'


            elif ((len(tokens) > 4) and (section_name == '#nonbond(9-6)')
                  and (pair_styles_selected & set(['class2','9-6','nonbond(9-6)']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_name = tokens[2]
                pair2ver[atom_name] = tokens[0]
                pair2ref[atom_name] = tokens[1]
                sigma = tokens[3]
                epsilon = tokens[4]
                pair2style[atom_name] = 'lj/class2/coul/long'
                pair_styles.add('lj/class2/coul/long')
                pair2params[atom_name] = (epsilon+' '+sigma)
                pair_mixing_style = 'sixthpower tail yes'
                #if pair_style_name.find('lj/class2') == 0:
                #    pair2params[atom_name] = (epsilon+' '+sigma)
                #    pair_mixing_style = 'sixthpower tail yes'


            elif (len(tokens) > 6) and (section_name == '#bond_increments'):
                if line.lstrip().find('!') == 0:
                    continue
                charge_pair_ver = tokens[0]
                charge_pair_ref = tokens[1]
                atom_names = map(EncodeAName, tokens[2:4])
                delta_q = [tokens[4], tokens[5]]
                if atom_names[0] > atom_names[1]:
                    atom_names.reverse()
                    delta_q.reverse()
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                charge_pair_priority[bond_name] = DeterminePriority(section_is_auto,
                                                                    tokens[2:4],
                                                                    float(charge_pair_ver[bond_name]))
                bond2chargepair[bond_name] = (delta_q[0] + ' ' + delta_q[1])


            elif ((len(tokens) > 5) and (section_name == '#quadratic_bond')
                  and (bond_styles_selected & set(['harmonic','quadratic','quadratic_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = DeterminePriority(section_is_auto,
                                                             tokens[2:4],
                                                             float(bond2ver[bond_name]))
                r0 = tokens[4]
                k = tokens[5]
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                bond2style[bond_name] = 'harmonic'
                bond2params[bond_name] = (k+' '+r0)


            elif ((len(tokens) > 6) and (section_name == '#morse_bond')
                  and (bond_styles_selected & set(['morse','morse_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = DeterminePriority(section_is_auto,
                                                             tokens[2:4],
                                                             float(bond2ver[bond_name]))
                r0 = tokens[4]
                D = tokens[5]
                alpha = tokens[6]
                sys.stderr.write('DEBUG: morse: atom_names = '+str(atom_names)+'\n')
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                bond2style[bond_name] = 'morse'
                bond2params[bond_name] = (D+' '+alpha+' '+r0)

            elif ((len(tokens) > 7) and (section_name == '#quartic_bond')
                  and (bond_styles_selected & set(['class2','quartic','quartic_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = DeterminePriority(section_is_auto,
                                                             tokens[2:4],
                                                             float(bond2ver[bond_name]))
                r0 = tokens[4]
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                K2 = tokens[5]
                K3 = tokens[6]
                K4 = tokens[7]
                bond2style[bond_name] = 'class2'
                bond2params[bond_name] = (r0+' '+K2+' '+K3+' '+K4)


        sys.stderr.write("parsing file pass3: look for (3-body) angle interactions...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')







            elif (len(tokens) > 6) and (section_name == '#quadratic_angle'):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:5]))
                angle_name = EncodeInteractionName(atom_names, section_is_auto)
                angle2ver[angle_name] = tokens[0]
                angle2ref[angle_name] = tokens[1]
                angle2priority[angle_name] = DeterminePriority(section_is_auto,
                                                               tokens[2:5],
                                                               float(angle2ver[angle_name]))
                theta0 = tokens[5]
                k = tokens[6]
                if not section_is_auto:
                    angle2theta0[angle_name] = theta0
                    sys.stderr.write('angle2theta0['+angle_name+'] = ' + str(theta0) + '\n')
                else:
                    angle2theta0_auto[(atom_names[0], atom_names[1], atom_names[2])] = theta0
                    sys.stderr.write('angle2theta0_auto['+str(atom_names)+'] = ' + str(theta0) + '\n')
                if (angle_styles_selected & set(['harmonic',
                                                 'quadratic',
                                                 'quadratic_angle'])):
                    angle2style[angle_name] = 'harmonic'
                    angle2params[angle_name] = (k+' '+theta0)
                elif (angle_styles_selected & set(['class2',
                                                   'quartic',
                                                   'quartic_angle'])):
                    # Then this is a special case of the class2 angle where
                    # the (theta-theta0)^3 and (theta-theta0)^4 terms = 0
                    angle2style[angle_name] = 'class2'
                    angle2params[angle_name] = (theta0+' '+k+' 0 0')



            elif ((len(tokens) > 8) and (section_name == '#quartic_angle')
                  and (angle_styles_selected & set(['class2','quartic','quartic_angle']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:5]))
                angle_name_sh = EncodeInteractionName(atom_names, section_is_auto)
                angle2ver[angle_name_sh] = tokens[0]
                angle2ref[angle_name_sh] = tokens[1]
                angle2priority[angle_name_sh] = \
                    DeterminePriority(section_is_auto,
                                      tokens[2:5],
                                      float(angle2ver[angle_name_sh]))
                theta0 = tokens[5]
                if not section_is_auto:
                    angle2theta0[angle_name_sh] = theta0
                    sys.stderr.write('angle2theta0['+angle_name_sh+'] = ' + str(theta0) + '\n')
                else:
                    angle2theta0_auto[(atom_names[0], atom_names[1], atom_names[2])] = theta0
                    sys.stderr.write('angle2theta0_auto['+str(atom_names)+'] = ' + str(theta0) + '\n')
                K2 = tokens[6]
                K3 = tokens[7]
                K4 = tokens[8]
                angle2style[angle_name_sh] = 'class2'
                angle2params_sh[angle_name_sh] = [theta0, K2, K3, K4]

            elif ((len(tokens) > 5) and
                  (section_name in ('#bond-bond', '#bond-angle')) and
                  (angle_styles_selected &
                   set(['class2', 'quartic', 'quartic_angle']))):
                if line.lstrip().find('!') == 0:
                    continue
                version = float(tokens[0])
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = map(EncodeAName, tokens[2:5])
                atom_names = ReverseIfEnds(aorig)
                K = ['', '']
                K[0] = tokens[5]
                K[1] = K[0]
                if len(tokens) > 6:
                    K[1] = tokens[6]
                order_reversed = aorig[0] > aorig[-1]
                if order_reversed:
                    K.reverse()
                if (section_name == '#bond-bond'):
                    angle2class2_bb_sh[angle_name] = K[0]
                elif (section_name == '#bond-angle'):
                    angle2class2_ba_sh[angle_name] = [k for k in K]









        sys.stderr.write("parsing file pass4: look for dihedrals(torsions) and impropers(out_of_plane)...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue


            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')




            elif (len(tokens) > 8) and (section_name == '#torsion_1'):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:6]))
                dihedral_name_sh = EncodeInteractionName(atom_names, section_is_auto)
                dihedral2ver_sh[dihedral_name_sh] = tokens[0]
                dihedral2ref_sh[dihedral_name_sh] = tokens[1]
                dihedral2priority_sh[dihedral_name_sh] = \
                    DeterminePriority(section_is_auto,
                                      tokens[2:6],
                                      float(dihedral2ver_sh[dihedral_name_sh]))
                K = tokens[6]
                n = tokens[7]
                d = tokens[8]
                w = '0.0'  #ignore: this is only used by the CHARMM force field

                if (dihedral_styles_selected & set(['charmm','torsion_1'])):
                    dihedral2style_sh[dihedral_name_sh] = 'charmm'
                    #dihedral2params_sh[dihedral_name_sh] = [K,n,d,w]
                    dihedral2params[dihedral_name_sh] = (K+' '+n+' '+d+' '+w)
                elif (dihedral_styles_selected & set(['class2','torsion_3'])):
                    # Then this is a special case of the class2 angle
                    # lacking the higher terms in the Fourier series
                    dihedral2style_sh[dihedral_name_sh] = 'class2'
                    dihedral2params_sh[dihedral_name_sh] = [K,d,0,0,0,0]
                      #= (K+' '+d+' '+
                      #    '0 0 '+'0 0')





            elif ((len(tokens) > 7) and (section_name == '#torsion_3')
                  and (dihedral_styles_selected & set(['class2','torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = ReverseIfEnds(map(EncodeAName, tokens[2:6]))
                dihedral_name_sh = EncodeInteractionName(atom_names, section_is_auto)
                dihedral2ver_sh[dihedral_name_sh] = tokens[0]
                dihedral2ref_sh[dihedral_name_sh] = tokens[1]
                dihedral2priority_sh[dihedral_name_sh] = \
                    DeterminePriority(section_is_auto,
                                      tokens[2:6],
                                      float(dihedral2ver_sh[dihedral_name_sh]))
                V1 = tokens[6]
                phi0_1 = tokens[7]
                V2 = phi0_2 = V3 = phi0_3 = '0.0'
                if len(tokens) > 9:
                    V2 = tokens[8]
                    phi0_2 = tokens[9]
                if len(tokens) > 11:
                    V3 = tokens[10]
                    phi0_3 = tokens[11]
                dihedral2style_sh[dihedral_name_sh] = 'class2'
                dihedral2params_sh[dihedral_name_sh] = [V1, phi0_1, V2, phi0_2, V3, phi0_3]




            elif ((len(tokens) > 6) and (section_name == '#middle_bond-torsion_3')
                  and (dihedral_styles_selected & set(['class2','torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                version = float(tokens[0])
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = map(EncodeAName, tokens[2:5])
                atom_names = ReverseIfEnds(aorig)

                Fmbt = [tokens[6], '0.0', '0.0']
                if len(tokens) > 7:
                    Fmbt[1] = tokens[7]
                if len(tokens) > 8:
                    Fmbt[2] = tokens[8]


                dihedral_name_sh = EncodeInteractionName(atom_names
                                                         section_is_auto)
                #sys.stderr.write('DEBUG: (a2,a3) = '+str((a2,a3))+', '
                #                 ' (b1,b2) = '+str(batoms)+'\n')
                #dihedral2ref_sh[dihedral_name_sh] = reference
                dihedral2style_sh[dihedral_name_sh] = 'class2'
                dihedral2class2_mbt_sh[dihedral_name_sh] = [F for F in Fmbt]
                dihedral2ver_mbt_sh[dihedral_name_sh] = version
                #dihedral2priority_sh[dihedral_name_sh] = DeterminePriority(section_is_auto,
                #                                                        tokens[2:6],
                #                                                        float(dihedral2ver[dihedral_name_sh]))




            elif ((len(tokens) > 6) and
                  (section_name in ('#end_bond-torsion_3',
                                    '#bond-bond_1_3')) and
                  (dihedral_styles_selected &
                   set(['class2', 'torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                version = float(tokens[0])
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = map(EncodeAName, tokens[2:6])
                atom_names = ReverseIfEnds(aorig)

                dihedral_name_sh = EncodeInteractionName(atom_names,
                                                         section_is_auto)

                dihedral2ref[dihedral_name_sh] = reference
                dihedral2style[dihedral_name_sh] = 'class2'
                if section_name == '#end_bond-torsion_3':
                    Febt = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
                    Febt[0][0] = tokens[6]
                    if len(tokens) > 7:
                        Febt[0][1] = tokens[7]
                    if len(tokens) > 8:
                        Febt[0][2] = tokens[8]
                    Febt[1][0] = Febt[0][0]
                    Febt[1][1] = Febt[0][1]
                    Febt[1][2] = Febt[0][2]
                    if len(tokens) > 9:
                        Febt[1][0] = tokens[9]
                    if len(tokens) > 10:
                        Febt[1][1] = tokens[10]
                    if len(tokens) > 11:
                        Febt[1][2] = tokens[11]
                    order_reversed = aorig[0] > aorig[-1]
                    if order_reversed:
                            Febt.reverse()
                    dihedral2class2_ebt_sh[dihedral_name_sh] = [ [F_ij for F_ij in F_i] for F_i in Febt] #deep copy of Febt[][]
                    dihedral2ver_ebt_sh[dihedral_name_sh] = version

                elif section_name == '#bond-bond_1_3':
                    Kbb13 = tokens[6]
                    #dihedral2ver_bb13[dihedral_name_sh] = version
                    dihedral2class2_bb13_sh[dihedral_name_sh] = Kbb13
                    dihedral2ver_bb13_sh[dihedral_name_sh] = version
                else:
                    assert(False)










            elif ((len(tokens) > 6) and
                  (section_name in ('#angle-torsion_3',
                                    '#angle-angle-torsion_1')) and
                  (dihedral_styles_selected &
                   set(['class2', 'torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                version = float(tokens[0])
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = map(EncodeAName, tokens[2:6])
                atom_names = ReverseIfEnds(aorig)

                dihedral_name_sh = EncodeInteractionName(atom_names,
                                                         section_is_auto)

                dihedral2ref[dihedral_name_sh] = reference
                dihedral2style[dihedral_name_sh] = 'class2'

                if section_name == '#angle-torsion_3':
                    Fat = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
                    Fat[0][0] = tokens[6]
                    if len(tokens) > 7:
                        Fat[0][1] = tokens[7]
                    if len(tokens) > 8:
                        Fat[0][2] = tokens[8]
                    Fat[1][0] = Fat[0][0]
                    Fat[1][1] = Fat[0][1]
                    Fat[1][2] = Fat[0][2]
                    if len(tokens) > 9:
                        Fat[1][0] = tokens[9]
                    if len(tokens) > 10:
                        Fat[1][1] = tokens[10]
                    if len(tokens) > 11:
                        Fat[1][2] = tokens[11]
                    order_reversed = aorig[0] > aorig[-1]
                    if order_reversed:
                        Fat.reverse()
                        Fat[0].reverse()
                        Fat[1].reverse()
                    dihedral2class2_at_sh = [ [F_ij for F_ij in F_i] for F_i in Fat] #deep copy of Fat[][]
                    dihedral2ver_at_sh[dihedral_name_sh] = version
                elif section_name == '#angle-angle-torsion_1':
                    Kaat = tokens[6]
                    dihedral2class2_aat_sh[dihedral_name_sh] = Kaat
                    dihedral2ver_aat_sh[dihedral_name_sh] = version
                else:
                    assert(False)












            elif ((len(tokens) > 8) and (section_name == '#out_of_plane')
                   and (improper_styles_selected & set(['cvff','out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names,_ignore  = OOPImproperNameSort(tokens[2:6])
                improper_name = EncodeInteractionName(atom_names, section_is_auto)
                improper2ver[improper_name] = tokens[0]
                improper2ref[improper_name] = tokens[1]
                improper2priority[improper_name] = DeterminePriority(section_is_auto,
                                                                     tokens[2:6],
                                                                     float(improper2ver[improper_name]))
                K = tokens[6]
                n = tokens[7]
                chi0 = tokens[8]
                improper2style[improper_name] = 'cvff'
                improper2params[improper_name] = (Kchi+' '+n+' '+chi0)
                improper_symmetry_subgraph = 'cenJswapIL'
                #if improper_style_name == 'cvff':
                #    improper2params[improper_name] = (Kchi+' '+n+' '+chi0)
                #    improper_symmetry_subgraph = 'cenJswapIL'

            elif ((len(tokens) > 7) and (section_name == '#wilson_out_of_plane')
                  and (improper_styles_selected & set(['class2','wilson_out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                sys.stderr.write('tokens = ' + str(tokens) + '\n')
                atom_names,_ignore = Class2ImproperNameSort(tokens[2:6])
                improper_name = EncodeInteractionName(atom_names, section_is_auto)
                improper2ver[improper_name] = tokens[0]
                improper2ref[improper_name] = tokens[1]
                improper2priority[improper_name] = DeterminePriority(section_is_auto,
                                                                     tokens[2:6],
                                                                     float(improper2ver[improper_name]))
                K = tokens[6]
                chi0 = tokens[7]
                improper2style[improper_name] = 'class2'
                #improper2class2_i[improper_name] = (K+' '+chi0)
                improper2params[improper_name] = (K+' '+chi0)
                improper_symmetry_subgraph = 'cenJsortIKL'
                #if improper_style_name == 'class2':
                #    improper2class2_i[improper_name] = (K+' '+chi0)
                #    improper_symmetry_subgraph = 'cenJsortIKL'

            elif ((len(tokens) > 6) and (section_name == '#angle-angle')
                  and (improper_styles_selected & set(['class2','wilson_out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names,_ignore = Class2ImproperNameSort(tokens[2:6])
                improper_name = EncodeInteractionName(atom_names, section_is_auto)
                improper2ver[improper_name] = tokens[0]
                improper2ref[improper_name] = tokens[1]
                improper2priority[improper_name] = DeterminePriority(section_is_auto,
                                                                     tokens[2:6],
                                                                     float(improper2ver[improper_name]))
                K = tokens[6]
                improper2cross[improper_name][ImCrossTermID(atom_names)] = K
                improper2style[improper_name] = 'class2'

            elif (len(tokens) > 0) and (section_name == '#out_of_plane-out_of_plane'):
                if line.lstrip().find('!') == 0:
                    continue
                display_OOP_OOP_warning = True

            elif (len(tokens) > 0) and (section_name == '#torsion-torsion_1'):
                if line.lstrip().find('!') == 0:
                    continue
                display_torsion_torsion_1_warning = True


            """
             --- these next few lines of code appear to be unnecessary.
             --- I'll probably delete this code in a later version
            elif (len(tokens) > 3) and (section_name == '#hbond_definition'):
                hbondID = tokens[1]
                if tokens[2] == 'distance':
                    hbond2distance[hbondID] = tokens[3]
                if tokens[2] == 'angle':
                    hbond2angle[hbondID] = tokens[3]
                if tokens[2] == 'donors':
                    hbond2donors[hbondID] = map(EncodeAName, tokens[2:])
                if tokens[2] == 'acceptors':
                    hbond2acceptors[hbondID] = map(EncodeAname(),tokens[2:])
            """


        if display_OOP_OOP_warning:
            sys.stderr.write('###########################################################\n'
                             'WARNING\n'
                             '      ALL \"out-of-plane_out-of_plane\" INTERACTIONS ARE IGNORED.\n'
                             '      CHECK THAT THESE TERMS ARE NEGLEGIBLY SMALL.\n'
                             '      \"out-of-plane_out-of_plane\" interactions are not yet supported in LAMMPS\n'
                             '      (...as of 2017-2-07)  There is no way that moltemplate can produce\n'
                             '      LAMMPS compatible parameter files for these interactions.\n'
                             '###########################################################\n')

        if display_torsion_torsion_1_warning:
            sys.stderr.write('###########################################################\n'
                             'WARNING\n'
                             '      ALL \"torsion_torsion_1\" INTERACTIONS ARE IGNORED.\n'
                             '      CHECK THAT THESE TERMS ARE NEGLEGIBLY SMALL.\n'
                             '      \"torsion_torsion_1\" interactions are not yet supported in LAMMPS\n'
                             '      (...as of 2017-2-07)  There is no way that moltemplate can produce\n'
                             '      LAMMPS compatible parameter files for these interactions.\n'
                             '###########################################################\n')


        sys.stderr.write(' done.\n'
                         'building lookup tables...')








        """
         --- these next few lines of code appear to be unnecessary.
         --- I'll probably delete them eventually
        if len(hbond2params) > 0:
            sys.stdout.write('\n\n  write_once("In Settings") {\n')
            if hbond_style == 'hbond/dreiding/lj':
                for hbondID, angle in hbond2angle:
                    hbond2params[hbondID] =  hbond2distance[hbondID]+' '+hbond2angle[hbondID]  ##<--this is not correct
            for hbondID, params in hbond2params:
                for donor in hbond2donors[hbondID]:
                    for acceptor in hbond2acceptors[hbondID]:
                        for hydrogen in hbond2hydrogens[hbondID]:
                            sys.stdout.write('pair_coeff @atom:'+donor+' @atom:'+acceptor+' '+hbond_style+' @atom:'+hydrogen+' i '+params+'\n')
            sys.stdout.write('  }   # (DREIDING style H-bond parameters)\n\n\n')
        """

        sys.stderr.write(" done.\n")
        sys.stderr.write("Trying all combinations of atom types...")





        ##################### POST-PROCESSING ########################





        for angle_name_sh in angle2class2_params_sh:
            assert(angle_name_sh in angle2class2_bb_sh)
            assert(angle_name_sh in angle2class2_ba_sh)

            is_auto = (angle_name_sh.find('auto_') == 0)

            atom_names = ExtractANames(angle_name_sh)

            num_angle_coeffs = 0

            atom_combos = [set([]), set([]), set([])]

            #*#atom_priorities = [{}, {}, {}]
            #*#atom_priorities[i][atom_name] = priority of i'th atom in interaction

            # We must consider every possible combination of atom types
            # which satisfy BOTH angle_equivalences and bond_equivalences.
            # ...AND we must consider BOTH regular AND auto equivalences.
            # For each combination generate a separate @angle interaction.
            # (I fear this will make the resulting .LT file large.)

            # Use different auto equivalence lookup tables for different
            # atoms in the interaction. (ie the "center" and "end" atoms)
            auto_angle2atom = [auto_angleend2atom,
                               auto_anglecenter2atom,
                               auto_angleend2atom]

            for i in range(0, 3):
                angle_atom_name = atom_names[i]
                sys.stderr.write('DEBUG: angle_atom_name = '+angle_atom_name+'\n')
                if not section_is_auto:
                    assert(angle_atom_name[-1] != '_')
                    # assume regular equivalences when looking up atom types
                    sys.stderr.write('DEBUG: equiv_angle2atom['+angle_atom_name+'] = '+
                                     str(equiv_angle2atom[angle_atom_name])+'\n')
                    for a in equiv_angle2atom[angle_atom_name]:
                        atom_combos[i].add(a)
                else:
                    assert((angle_atom_name[-1] == '_') or (ange_atom_name[0] == '*'))
                    # assume "auto" equivalences when looking up atom types
                    sys.stderr.write('DEBUG: auto_angle2atom['+str(i)+']['+angle_atom_name+'] = \n'
                                     '       '+str(equiv_angle2atom[i][angle_atom_name])+'\n')
                    for a in auto_angle2atom[i][angle_atom_name]:
                        atom_combos[i].add(a)

            found_at_least_one = False
            #*#for a1, a1priority in atom_priorities[0].items():
            #*#    for a2, a2priority in atom_priorities[1].items():
            #*#        for a3, a3priority in atom_priorities[2].items():
            for a1 in atom_combos[0]:
                for a2 in atom_combos[1]:
                    for a3 in atom_combos[2]:
                        #sys.stderr.write('atom2auto_bond = '+str(atom2auto_bond)+'\n')
                        bond_data1 = LookupBondLength(a1, a2,
                                                      atom2equiv_bond,
                                                      bond2r0,
                                                      atom2auto_bond,
                                                      bond2r0_auto)
                        if bond_data1 != None:
                            # Save time by only continuing if a bond was
                            # found between a1 and a2
                            bond_data2 = LookupBondLength(a2, a3,
                                                          atom2equiv_bond,
                                                          bond2r0,
                                                          atom2auto_bond,
                                                          bond2r0_auto)
                        if (bond_data1 and bond_data2):
                            #bond lengths:
                            r0s = [0.0, 0.0]
                            #equivalent atom names used to lookup the bonds:
                            batoms = [['', ''], ['', '']]
                            r0s[0], batoms[0] = bond_data1
                            r0s[1], batoms[1] = bond_data2
                            found_at_least_one = True
                            angle_name_full = EncodeInteractionName(atom_names +
                                                                    batoms[0] + batoms[1],
                                                                    #+ [str(r0s[0]),
                                                                    #   str(r0s[1])],
                                                                    section_is_auto)
                            #sys.stderr.write('DEBUG: (a1,a2,a3) = '+str((a1,a2,a3))+', '
                            #                 ' (b11,b12,b21,b22) = '+str(batoms)+'\n')
                            angle2ref[angle_name_full] = reference
                            angle2style[angle_name_full] = 'class2'
                            Kbb = angle2class2_bb_sh[angle_name_sh]
                            Kba = angle2class2_ba_sh[angle_name_sh]
                            angle2class2_bb[angle_name_full] = (Kbb+' '+r0s[0]+' '+r0s[1])
                            angle2priority_bb = DeterminePriority(is_auto,
                                                                  batoms[0] + batoms[1],
                                                                  angle2class2ver_bb_sh[angle_name_sh])
                            angle2class2_ba[angle_name_full] = (Kba[0]+' '+Kba[1]+' '+r0s[0]+' '+r0s[1])
                            angle2sym_ba = (r0s[0] == r0s[1])
                            angle2priority_ba = DeterminePriority(is_auto,
                                                                  batoms[0] + batoms[1],
                                                                  angle2class2ver_ba_sh[angle_name_sh])
                            version = max((angle2ver_sh[angle_name_sh],
                                           angle2class2ver_bb_sh[angle_name_sh],
                                           angle2class2ver_ba_sh[angle_name_sh]))
                            angle2ver[angle_name_full] = version
                            #angle2priority[angle_name] = \
                            #    DeterminePriority(section_is_auto,
                            #                      aorig,
                            #                      float(angle2ver[angle_name]))
                            angle2priority[angle_name_full] = \
                               max((angle2class2priority[angle_name_sh],
                                    angle2class2priority_bb,
                                    angle2class2priority_ba))

                            if num_angle_coeffs < len(angle2class2_bb):
                                sys.stderr.write('DEBUG: '+section_name[1:]+' r0 ('+angle_name+') = ('+r0s[0]+', '+r0s[1]+')\n')
                                sys.stderr.write('DEBUG: num_angle_coeffs_bb = len(angle2class2_bb) = '+str(len(angle2class2_bb))+'\n')
                                sys.stderr.write('DEBUG: '+section_name[1:]+' r0 ('+angle_name+') = ('+r0s[0]+', '+r0s[1]+')\n')
                                sys.stderr.write('DEBUG: num_angle_coeffs_ba = len(angle2class2_ba) = '+str(len(angle2class2_ba))+'\n')
                            num_angle_coeffs = len(angle2class2_bb)

                        if ((not angle2sym_ba)
                            and
                            (anames[0] == anames[3])):
                            raise InputError('Error: Unsupported dihedral interaction: \"@angle:'+str(angle_name_sh)+'\"\n'
                                             '       This interaction has symmetric atom names:\n'
                                             ', '.join(anames)+'\n'
                                             '       and yet it lacks symmetry in the corresponding force field parameters.\n'
                                             '       (If this is not a mistake in the .frc file, then explain\n'
                                             '        why to andrew so he can fix this.)\n')


            if not found_at_least_one:
                #raise InputError('Error: Undefined bonds for bond-bond interactions:\n'
                #                 '       '+str(atom_names)+'\n')
                sys.stderr.write('WARNING: Undefied bond length for ' +
                                 #'         '+
                                 section_name[1:] + ' interaction: ' +
                                 ' '.join(atom_names)+'\n')
            #sys.stderr.write('bond_names = ' + str(bond_names) + '\n')





       ############ POST-PROCESSING DIHEDRALS ###########



        for dihedral_name_sh in dihedral2class2_params_sh:
            assert(dihedral_name_sh in dihedral2class2_mbt_sh)
            assert(dihedral_name_sh in dihedral2class2_ebt_sh)
            assert(dihedral_name_sh in dihedral2class2_bb13_sh)
            assert(dihedral_name_sh in dihedral2class2_at_sh)
            assert(dihedral_name_sh in dihedral2class2_aat_sh)

            is_auto = (dihedral_name_sh.find('auto_') == 0)

            atom_names = ExtractANames(dihedral_name_sh) 

            num_dihedral2class2 = 0

            atom_combos = [set([]), set([]), set([]), set([])]

            #*#atom_priorities = [{}, {}, {}, {}]
            #*#atom_priorities[i][atom_name] = priority of i'th atom in interaction

            # We must consider every possible combination of atom types
            # which satisfy all three:
            #   dihedral_equivalences
            #   bond_equivalences
            #   angle_equivalences
            # ...AND we must consider BOTH regular AND auto equivalences.
            # For each combination generate a separate @dihedral interaction.
            # (I fear this will make the resulting .LT file large.)

            # Use different auto equivalence lookup tables for different
            # atoms in the interaction. (ie the "center" and "end" atoms)
            auto_dihedral2atom = [auto_dihedralend2atom,
                                  auto_dihedralcenter2atom,
                                  auto_dihedralcenter2atom,
                                  auto_dihedralend2atom]

            for i in range(0, 4):
                dihedral_atom_name = atom_names[i]
                sys.stderr.write('DEBUG: dihedral_atom_name = '+dihedral_atom_name+'\n')
                if not section_is_auto:
                    assert(dihedral_atom_name[-1] != '_')
                    # assume regular equivalences when looking up atom types
                    sys.stderr.write('DEBUG: equiv_dihedral2atom['+dihedral_atom_name+'] = '+
                                     str(equiv_dihedral2atom[dihedral_atom_name])+'\n')
                    for a in equiv_dihedral2atom[dihedral_atom_name]:
                        atom_combos[i].add(a)
                else:
                    assert((dihedral_atom_name[-1] == '_') or (ange_atom_name[0] == '*'))
                    # assume "auto" equivalences when looking up atom types
                    sys.stderr.write('DEBUG: auto_dihedral2atom['+str(i)+']['+dihedral_atom_name+'] = \n'
                                     '       '+str(equiv_dihedral2atom[i][dihedral_atom_name])+'\n')
                    for a in auto_dihedral2atom[i][dihedral_atom_name]:
                        atom_combos[i].add(a)

            found_at_least_one = False
            #*#for a1, a1priority in atom_priorities[0].items():
            #*#    for a2, a2priority in atom_priorities[1].items():
            #*#        for a3, a3priority in atom_priorities[2].items():
            #*#            for a4, a3priority in atom_priorities[3].items():
            for a1 in atom_combos[0]:
                for a2 in atom_combos[1]:
                    #sys.stderr.write('atom2auto_bond = '+str(atom2auto_bond)+'\n')
                    bond_data12 = LookupBondLength(a1, a2,
                                                   atom2equiv_bond,
                                                   bond2r0,
                                                   atom2auto_bond,
                                                   bond2r0_auto)
                    if bond_data12 == None:
                        # Save time by only continuing if a bond was
                        # found between a1 and a2
                        continue
                    for a3 in atom_combos[2]:
                        bond_data23 = LookupBondLength(a2, a3,
                                                       atom2equiv_bond,
                                                       bond2r0,
                                                       atom2auto_bond,
                                                       bond2r0_auto)
                        if bond_data23 == None:
                            # Save time by only continuing if a bond was
                            # found between a2 and a3
                            continue

                        angle_data123 = LookupRestAngle(a1, a2, a3,
                                                        atom2equiv_angle,
                                                        angle2theta0,
                                                        [atom2auto_angleend,
                                                         atom2auto_anglecenter,
                                                         atom2auto_anglecenter],
                                                        angle2theta123_auto)
                        if angle_data123 != None:
                            # Save time by only continuing if a angle was
                            # found between a1, a2, a3
                            continue

                        for a4 in atom_combos[3]:
                            bond_data34 = LookupBondLength(a3, a4,
                                                           atom2equiv_bond,
                                                           bond2r0,
                                                           atom2auto_bond,
                                                           bond2r0_auto)
                            if bond_data34 == None:
                                # Save time by only continuing if a bond was
                                # found between a3 and a4
                                continue

                            #rest bond lengths:
                            r0s = [0.0, 0.0, 0,0]
                            #equivalent atom names used to lookup the bonds:
                            batoms = [['', ''], ['', ''], ['','']]
                            r0s[0], batoms[0] = bond_data12
                            r0s[1], batoms[1] = bond_data23
                            r0s[2], batoms[2] = bond_data34

                            angle_data234 = LookupRestAngle(a2, a3, a4,
                                                            atom2equiv_angle,
                                                            angle2theta0,
                                                            [atom2auto_angleend,
                                                             atom2auto_anglecenter,
                                                             atom2auto_anglecenter],
                                                            angle2theta234_auto)
                            if angle_data234 != None:
                                # Save time by only continuing if a angle was
                                # found between a2, a3, a4
                                continue

                            #rest angles:
                            theta0s = [0.0, 0.0]
                            #equivalent atom names used to lookup angles:
                            aatoms = [['', '',''], ['', '','']]
                            theta0s[0], aatoms[0] = angle_data123
                            theta0s[1], aatoms[1] = angle_data234
                            found_at_least_one = True
                            order_reversed = aorig[0] > aorig[-1]
                            if order_reversed:
                                theta0s.reverse()
                                aatoms.reverse()
                                aatoms[0].reverse()
                                aatoms[1].reverse()

                            dihedral_name_full = dihedral_name_sh + \
                                                 EncodeInteractionName(batoms[0] + batoms[1] + batoms[2] +
                                                                       aatoms[0] + aatoms[1],
                                                                       False)
                            else:
                                assert(batoms[0][1] == batoms[1][0])
                                assert(batoms[1][1] == batoms[2][0])
                                assert(aatoms[0][1] == aatoms[1][0])
                                assert(aatoms[0][2] == aatoms[1][1])
                                dihedral_name_full = dihedral_name_sh + 
                                                     EncodeInteractionName([batoms[0][0], batoms[0][1]
                                                                            batoms[2][0], batoms[2][1],
                                                                            aatoms[0][0], aatoms[0][1],
                                                                            aatoms[0][2], aatoms[1][0]],
                                                                           False)

                            found_at_least_one = True

                            ########### Fourier terms ###########
                            V_phi0_params = dihedral2params_sh[dihedral_name_sh]
                            dihedral2params_sh[dihedral_name_full] = ' '.join(V_phi0_params)

                            ########### "mbt", "ebt", and "aat" terms ###########
                            # "mbt" terms:
                            Fmbt = dihedral2class2_mbt_sh[dihedral_name_sh]
                            dihedral2class2_mbt[dihedral_name_full] = \
                               (Fmbt[0]+' '+Fmbt[1]+' '+Fmbt[2]+' '+r0s[1])
                            dihedral2priority_mbt = DeterminePriority(is_auto,
                                                                      batoms[1],
                                                                      dihedral2class2ver_mbt[dihedral_name_sh])

                            # "ebt" terms:
                            Febt = dihedral2class2_ebt_sh[dihedral_name_sh]
                            dihedral2class2_ebt[dihedral_name_full]= (Febt[0][0] + ' ' +
                                                                      Febt[0][1] + ' ' +
                                                                      Febt[0][2] + ' ' +
                                                                      Febt[1][0] + ' ' +
                                                                      Febt[1][1] + ' ' +
                                                                      Febt[1][2] + ' ' +
                                                                      r0s[0]+' '+r0s[2])
                            dihedral2sym_ebt = ((Febt[0][0] == Febt[1][0]) and
                                                (Febt[0][1] == Febt[1][1]) and
                                                (Febt[0][2] == Febt[1][2]) and
                                                (r0s[0] == r0s[2]))
                            dihedral2priority_ebt = DeterminePriority(is_auto,
                                                                      batoms[0] + batoms[2],
                                                                      dihedral2class2ver_ebt[dihedral_name_sh])

                            #(Note:  large atom_priority number <==> low priority
                            # Only one of the atom priority numbers should be > 0)

                            # "bb13" terms:
                            Kbb13 = dihedral2class2_bb13_sh[dihedral_name_sh]
                            dihedral2class2_bb13[dihedral_name_full] = (Kbb13+' '+r0s[0]+' '+r0s[2])
                            dihedral2sym_bb13 = (r0s[0] == r0s[2])
                            dihedral2priority_bb13 = DeterminePriority(is_auto,
                                                                       batoms[0] + batoms[2],
                                                                       dihedral2class2ver_bb13[dihedral_name_sh])


                            ########### "at" and "aat" terms ###########
                            # "at" terms:
                            Fat = dihedral2class2_at_sh[dihedral_name_sh]
                            dihedral2class2_at[dihedral_name_full] = \
                                (Fat[0][0] + ' ' +
                                 Fat[0][1] + ' ' +
                                 Fat[0][2] + ' ' +
                                 Fat[1][0] + ' ' +
                                 Fat[1][1] + ' ' +
                                 Fat[1][2] + ' ' +
                                 theta0s[0] + ' ' +
                                 theta0s[1])
                            dihedral2sym_at = ((Fat[0][0] == Fat[1][0]) and
                                               (Fat[0][1] == Fat[1][1]) and
                                               (Fat[0][2] == Fat[1][2]) and
                                               (theta0[0] == theta0[1]))
                            dihedral2priority_at = DeterminePriority(is_auto,
                                                                     aatoms[0] + aatoms[1],
                                                                     dihedral2class2ver_at[dihedral_name_sh])
                            # "aat" terms:
                            Kaat = dihedral2class2_aat_sh[dihedral_name_sh]
                            dihedral2class2_aat[dihedral_name_full] = \
                                               (Kaat+' '+
                                                theta0s[0]+' '+
                                                theta0s[1])
                            dihedral2sym_aat[dihedral_name_full] = (theta0[0] == theta0[1])
                            dihedral2priority_aat = DeterminePriority(is_auto,
                                                                      aatoms[0] + aatoms[1],
                                                                      dihedral2class2ver_aat[dihedral_name_sh])

                            if len(dihedral2class2_ebt) > num_dihedral2class2:
                                sys.stderr.write('DEBUG: dihedral['+dihedral_name_full+']:\n'
                                                 '(r12,r23,r34) = ('
                                                 +r0s[0]+','+r0s[1]+','+r0s[2]+') \n'
                                                 '(theta123,theta234) = ('
                                                 +theta0s[0]+','+theta0s[1]+') \n')
                                sys.stderr.write('DEBUG: num_dihedral_coeffs = len(dihedral2class2_ebt) = '
                                                 +str(len(dihedral2class2_ebt))+'\n')
                        version = max((dihedral2ver_sh[dihedral_name_sh],
                                       dihedral2class2ver_mbt_sh[dihedral_name_sh],
                                       dihedral2class2ver_ebt_sh[dihedral_name_sh],
                                       dihedral2class2ver_bb13_sh[dihedral_name_sh],
                                       dihedral2class2ver_at_sh[dihedral_name_sh],
                                       dihedral2class2ver_aat_sh[dihedral_name_sh]))

                        dihedral2class2ver[dihedral_name_full] = version
                        dihedral2priority[dihedral_name_full] = \
                            max((dihedral2class2priority[dihedral_name_sh],
                                 dihedral2class2priority_mbt,
                                 dihedral2class2priority_ebt,
                                 dihedral2class2priority_bb13,
                                 dihedral2class2priority_at,
                                 dihedral2class2priority_aat)

                        num_dihedral2class2 = len(dihedral2class2_ebt)

                        if ((not (dihedral2sym_ebt and
                                  #dihedral2sym_mbt and
                                  # (note: symmetry doesn't make sense for mbt)
                                  dihedral2sym_at and
                                  dihedral2sym_aat and
                                  dihedral2sym_bb13))
                            and
                            ((anames[0] == anames[3]) and
                             (anames[1] == anames[2]))):
                            raise InputError('Error: Unsupported dihedral interaction: \"@dihedral:'+str(dihedral_name_sh)+'\"\n'
                                             '       This interaction has symmetric atom names:\n'
                                             ', '.join(anames)+'\n'
                                             '       and yet it lacks symmetry in the corresponding force field parameters.\n'
                                             '       (If this is not a mistake in the .frc file, then explain\n'
                                             '        why to andrew so he can fix this.)\n')



            #sys.stderr.write('DEBUG: number of interactions = '+str(len(dihedral2class2_bb))+'\n')
            if not found_at_least_one:
                #raise InputError('Error: Undefined bonds for bond-bond interactions:\n'
                #                 '       '+str(atom_names)+'\n')
                sys.stderr.write('WARNING: Undefined bond length (r0) or rest angle (theta0) for\n'+
                                 #'         '+
                                 + ' dihedral interaction between these atoms: ' +
                                 ' '.join(atom_names)+'\n')
                #sys.stderr.write('bond_names = ' + str(bond_names) + '\n')









       ############ POST-PROCESSING IMPROPERS ###########




        # Collect information from the different terms in a class2 improper:
        # http://lammps.sandia.gov/doc/improper_class2.html

        for improper_name in improper2cross:
            # Loop over the neighbors of the central atom in each improper
            # interaction and collect all the Mi and Ti parameters. Collect 
            # them in the order they appear in the formula for the Eaa
            # term as it appears in the documentation for improper_style class2:
            # 
            #    http://lammps.sandia.gov/doc/improper_class2.html
            #
            # Eaa = M1 (Tijk - T0)(Tkjl - T2) +   #common leaf node: k (index 2)
            #       M2 (Tijk - T0)(Tijl - T1) +   #common leaf node: i (index 0)
            #       M3 (Tijl - T1)(Tkjl - T2)     #common leaf node: l (index 3)
            # (I'm trying to match the variable names used in this web page
            #  I wish the author had chosen the M1,M2,M3, T1,T2,T3 order in more
            #  symmetric way, or at least in a way that makes more sense to me.)

            is_auto = IsAutoInteraction(improper_name)   # is this an "auto" interaction?

            anames = ExtractANames(improper_name)        # names of all 4 atoms
            lnames = [anames[0], anames[2], anames[3]]   # names of "leaf" atoms

            #M1     = improper2cross[improper_name][ 2 ]
            #M2     = improper2cross[improper_name][ 0 ]
            #M3     = improper2cross[improper_name][ 3 ]
            M1     = improper2cross[improper_name][ImCrossTermID([anames[0],
                                                                  anames[1],
                                                                  anames[2],
                                                                  anames[3]])]
            M2     = improper2cross[improper_name][ImCrossTermID([anames[2],
                                                                  anames[1],
                                                                  anames[0],
                                                                  anames[3]])]
            M3     = improper2cross[improper_name][ImCrossTermID([anames[0],
                                                                  anames[1],
                                                                  anames[3],
                                                                  anames[2]])]

            angle_name_l = ReverseIfEnds([anames[0], anames[1], anames[2]])
            angle_name = EncodeInteractionName(angle_name_l, is_auto)
            theta01 = angle2theta0[angle_name]

            angle_name_l = ReverseIfEnds([anames[0], anames[1], anames[3]])
            angle_name = EncodeInteractionName(angle_name_l, is_auto)
            theta02 = angle2theta0[angle_name]

            angle_name_l = ReverseIfEnds([anames[2], anames[1], anames[3]])
            angle_name = EncodeInteractionName(angle_name_l, is_auto)
            theta03 = angle2theta0[angle_name]

            improper2class2_aa[improper_name] = (M1 + ' ' + M2 + ' ' + M3 + ' '+
                                                 theta01 + ' ' +
                                                 theta02 + ' ' +
                                                 theta03)

            # ###### Symmetry: ######
            # Unfortunately, it's time to wade into the messy issue of symmetry.
            #    We desire a way to detect whether an improper interaction
            # between 4 atoms is invariant with respect to atom reordering
            # of the 3 peripheral "leaf" atoms which surround the central atom.
            # In principle, any rearrangement of atoms would require a separate
            # class2 improper interaction.  However, in some cases, when the
            # parameters for these rearrangements are symmetric, we can detect
            # that and warn moltemplate that it is not necessary to generate new
            # improper interactions for every conceivable permutation of these
            # atoms.  Figuring out when it is safe to do that is a headache.
            #   (...but it's necessary.  Otherwise each junction in the molecule
            #   will generate 3*2*1=6 improper interactions which are usually
            #   redundant.  This will slow down the simulation significantly
            #   and may make it difficult to compare the resulting LAMMPS 
            #   input files with those generated by other tools like msi2lmp.)
            #
            # To make this easier, I store the parameters in arrays which 
            # are arranged in a more symmetric way
            M = [0.0, 0.0, 0.0]
            theta0 = [0.0, 0.0, 0.0]
            # noti3[i] = the sorted tuple of integers from the 
            #            set {0,1,2} which remain after deleting i
            noti3 = ((1,2), (0,2), (0,1))
            i_neigh = [ ([0,2,3][ noti3[i][0] ],   # neighbor leaves of ith leaf
                         [0,2,3][ noti3[i][1] ]) for i in range(0,3)]
            for i in range(0, 3):
                # You will notice the pattern "[0,2,3][i]" appears often in the
                # code below because for class 2 force-fields, the second atom
                # (with index 1) is the central atom ("hub" atom), and the three
                # that surround it ("leaf" atoms) have indices 0,2,3.  I want
                # to skip over the central atoms and loop over the leaf atoms
                imTermID = ImCrossTermID([anames[ i_neigh[i][0] ],
                                          anames[ 1 ],
                                          anames[ [0,2,3][i] ],
                                          anames[ i_neigh[i][1] ]])
                M[i] = float(improper2cross[improper_name][imTermID])
                #i_leaf = [0,2,3][i]
                #M[i] = float(improper2cross[improper_name][ i_leaf ])
                angle_name_l = ReverseIfEnds([anames[i_neigh[i][0]],
                                              anames[i],
                                              anames[i_neigh[i][1]]])
                angle_name = EncodeInteractionName(angle_name_l, is_auto)
                theta0[i] = float(angle2theta0[angle_name])

            for i in range(0, 3):
                if ((M[i_neigh[i][0]] == M[i_neigh[i][1]]) and
                    (theta0[ i_neigh[i][1] ] == theta0[ i_neigh[i][1] ])):
                    # Then it is safe to swap the order of these two atoms in
                    # the list of atoms when looking up force-field parameters
                    improper2sym[improper_name].add(i_neigh[i][0])
                    improper2sym[improper_name].add(i_neigh[i][1])
                    # Later, I can use these to decide whether or not I need to
                    # change the default script with symmetry rules. (I'm hoping
                    # that "cenJsortIKL.py" should work in most cases.)
                else:
                    if anames[i_neigh[i][0]] == anames[i_neigh[i][1]]:
                        raise InputError('Error: Unsupported improper interaction: \"@improper:'+str(improper_name)+'\"\n'
                                         '       This interaction has matching aton aliases:\n'
                                         '       (@atom:'+str(anames[i_neigh[i][0]])+
                                         ', @atom:'+str(anames[i_neigh[i][1]])+')\n'
                                         '       and yet it lacks symmetry in the corresponding force field parameters.\n'
                                         '       (If this is not a mistake in the .frc file, then ask andrew to\n'
                                         '       fix this limitation.)\n')









        sys.stderr.write("done\n")
        sys.stderr.write("Converting to moltemplate format...\n")





        ##################### BEGIN WRITING FILE #####################





        sys.stdout.write("# This file was generated automatically using:\n")
        sys.stdout.write("# " + g_program_name + " " + " ".join(sys.argv[1:]) + "\n")
        sys.stdout.write("\n\n")
        sys.stdout.write(ffname + " {\n\n")
        
        sys.stdout.write("  # Below we will use lammps \"set\" command to assign atom charges\n"
                         "  # by atom type.  http://lammps.sandia.gov/doc/set.html\n\n")
        
        sys.stdout.write("\n"
                         "  #        AtomType    Mass     # \"Description\" (version, reference)\n\n")
        sys.stdout.write("  write_once(\"Data Masses\") {\n")
        for atype in atom2mass:
            sys.stdout.write("    @atom:" + atype + " " + str(atom2mass[atype]))
            sys.stdout.write("  # ")
            if atype in atom2elem:
                sys.stdout.write(atom2elem[atype] + ", ")
            #sys.stdout.write(atom2descr[atype])
            sys.stdout.write("\"" + atom2descr[atype] + "\"")
            sys.stdout.write(" (")
            if atype in atom2numbonds:
                sys.stdout.write("nbonds="+str(atom2numbonds[atype])+", ")
            sys.stdout.write("ver=" + atom2ver[atype] +
                             ", ref=" + atom2ref[atype])
            sys.stdout.write(")\n")
        sys.stdout.write("  } #(end of atom masses)\n\n\n")












        sys.stdout.write("  # ---------- EQUIVALENCE CATEGORIES for bonded interaction lookup ----------\n"
                         "  #   Each type of atom has a separate ID used for looking up bond parameters\n"
                         "  #   and a separate ID for looking up 3-body angle interaction parameters\n"
                         "  #   and a separate ID for looking up 4-body dihedral interaction parameters\n"
                         "  #   and a separate ID for looking up 4-body improper interaction parameters\n"
                         #"  #   (This is because there are several different types of sp3 carbon atoms\n"
                         #"  #   which have the same torsional properties when within an alkane molecule,\n"
                         #"  #   for example.  If they share the same dihedral-ID, then this frees us\n"
                         #"  #   from being forced define separate dihedral interaction parameters\n"
                         #"  #   for all of them.)\n"
                         "  #   The complete @atom type name includes ALL of these ID numbers.  There's\n"
                         "  #   no need to force the end-user to type the complete name of each atom.\n"
                         "  #   The \"replace\" command used below informs moltemplate that the short\n"
                         "  #   @atom names we have been using abovee are equivalent to the complete\n"
                         "  #   @atom names used below:\n\n")

        for atype in atom2ffid:
            #ffid = atype + "_ffid" + atom2ffid[atype]
            sys.stdout.write("  replace{ @atom:" + atype +
                             " @atom:" + atom2ffid[atype] + " }\n")
        
        sys.stdout.write("\n\n\n\n")
        
        
        sys.stdout.write("  # --------------- Non-Bonded Interactions: ---------------------\n"
                         "  # Syntax:\n"
                         "  # pair_coeff    AtomType1    AtomType2   pair_style_name  parameters...\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
                        
        for atype in pair2params:
            assert(atype in pair2style)
            assert(atype in atom2equiv_pair)
        
            sys.stdout.write("    pair_coeff " +
                             "@atom:" + atom2equiv_pair[atype] + " " + 
                             "@atom:" + atom2equiv_pair[atype] + " " + 
                             pair2style[atype] + " " +
                             pair2params[atype] +
                             " # (ver=" + pair2ver[atype] +
                             ", ref=" +pair2ref[atype] + ")\n")
        sys.stdout.write("  } #(end of pair_coeffs)\n\n\n\n")
        






        ################# Print Charge By Bond Interactions ##################
        charge_pair_priority_high_to_low = sorted(charge_pair_priority.items(),
                                                  key=itemgetter(1),
                                                  reverse=True)

        if len(charge_pair_priority_high_to_low) > 0:
            sys.stdout.write("  # ---------- Charge By Bond (a.k.a. \"bond equivalences\") ----------\n")
            sys.stdout.write("  write_once(\"Data Charge By Bond\") {\n")
            # Print rules for generating (2-body) "bond" interactions:
            sys.stdout.write('\n\n\n'
                             '  write_once("Charge By Bond") {\n')
            for bond_name, in charge_pair_priority_high_to_low:
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(bond_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(bond_name):
                    if include_auto_equivalences:
                        sys.stdout.write(' @atom:*,aq' + anames[0] +
                                         ',ab*,aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' @atom:*,aq' + anames[1] +
                                         ',ab*,aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' ' + bond2chargepair[bond_name] +
                                         " # (ver=" + charge_pair_ver[bond_name] +
                                         ", ref=" + charge_pair_ref[bond_name] + ")\n")
                    else:
                        continue
                else:
                    sys.stdout.write(' @atom:' + anames[0] + 'b*,a*,d*,i* ' +
                                     ' @atom:' + anames[1] + 'b*,a*,d*,i* ' +
                                     ' ' + bond2chargepair[bond_name] +
                                     " # (ver=" + charge_pair_ver[bond_name] +
                                     ", ref=" + charge_pair_ref[bond_name] + ")\n")
                    #  --- or should I use this instead ? ----
                    # (I'm not sure how bond_increments use equivalences)
                    #sys.stdout.write(' @atom:*,b' + anames[0] + ',a*,d*,i* ' +
                    #                 ' @atom:*,b' + anames[1] + ',a*,d*,i* ' +
                    #                 ' ' + bond2chargepair[bond_name] +
                    #                 '\n')
            sys.stdout.write('  } #(end of Charge by Bond (bond equivalences))\n\n'
                             '\n\n\n\n')







        ################# Print 2-body Bond Interactions ##################

        bond_names_priority_high_to_low = sorted(bond2priority.items(),
                                                 key=itemgetter(1),
                                                 reverse=True)

        if len(bond_names_priority_high_to_low) > 0:
            sys.stdout.write("  # --------------- Bond Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (2-body) "bond" interactions: --\n'
                             '  #  BondType  AtomType1  AtomType2\n')
            sys.stdout.write('\n'
                             '  write_once("Data Bonds By Type") {\n')
            for bond_name, in bond_names_priority_high_to_low:
                if not (bond2style[bond_name] in
                        bond_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(bond_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(bond_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @bond:' + bond_name + ' ' +
                                         ' @atom:*,aq*,ab' + anames[0] +
                                         ',aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' @atom:*,aq*,ab' + anames[1] +
                                         ',aae*,aac*,ade*,adc*,aie*,aic*' +
                                         '\n')
                    else:
                        continue
                else:
                    sys.stdout.write('    @bond:' + bond_name + ' ' +
                                     ' @atom:*,b' + anames[0] + ',a*,d*,i* ' +
                                     ' @atom:*,b' + anames[1] + ',a*,d*,i* ' +
                                     '\n')

            sys.stdout.write('  }  # end of "Data Bonds By Type" section\n'
                             '\n')

            # Print the force-field parameters for these bond interactions:
            sys.stdout.write('\n\n'
                             '  # ------------ Bond Parameters: ----------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write('    # '+bond_style2docs[bond_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # bond_coeff BondTypeName  BondStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for bond_name, in bond_names_priority_high_to_low:
                if not (bond2style[bond_name] in
                        bond_styles_selected):
                    continue
                # Did the user ask us to include "auto" interactions?
                if (IsAutoInteraction(bond_name) and
                    (not include_auto_equivalences)):
                    continue
                sys.stdout.write('    bond_coeff @bond:'+bond_name+' '+
                                 bond2style[bond_name] +
                                 bond2params[bond_name] +
                                 " # (ver=" + bond2ver[bond_name] +
                                 ", ref=" +bond2ref[bond_name] + ")\n")

            sys.stdout.write('  }  # end of bond_coeff commands\n'
                             '\n\n')






        ################# Print 3-body Angle Interactions ##################

        angle_names_priority_high_to_low = sorted(angle2priority.items(),
                                                  key=itemgetter(1),
                                                  reverse=True)

        if len(angle_names_priority_high_to_low) > 0:
            sys.stdout.write("  # --------------- Angle Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (3-body) "angle" interactions: --\n'
                             '  #  AngleType AtomType1 AtomType2 AtomType3  [BondType1 BondType2]\n')
            sys.stdout.write('\n'
                             '  write_once("Data Angles By Type") {\n')
            for angle_name, in angle_names_priority_high_to_low:
                if not (angle2style[angle_name] in
                        angle_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(angle_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(angle_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @angle:' + angle_name + ' ' +
                                         ' @atom:*,aq*,ab*,aae' + anames[0] +
                                         ',aac*,ade*,adc*,aie*,aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac'+anames[1] +
                                         ',ade*,adc*,aie*,aic*' +
                                         ' @atom:*,aq*,ab*,aae' + anames[2] +
                                         ',aac*,ade*,adc*,aie*,aic*' +
                                         '\n')
                    else:
                        continue
                else:
                    sys.stdout.write('    @angle:' + angle_name + ' ' +
                                     ' @atom:*,b*,a' + anames[0] + ',d*,i* ' +
                                     ' @atom:*,b*,a' + anames[1] + ',d*,i* ' +
                                     '\n')

            sys.stdout.write('  }  # end of "Data Angles By Type" section\n'
                             '\n')

            # Print the force-field parameters for these angle interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Angle Force Field Parameters: -------')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write('    # '+angle_style2docs[angle_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # angle_coeff AngleTypeName  AngleStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for angle_name, in angle_names_priority_high_to_low:
                if not (angle2style[angle_name] in
                        angle_styles_selected):
                    continue
                # Did the user ask us to include "auto" interactions?
                if (IsAutoInteraction(angle_name) and
                    (not include_auto_equivalences)):
                    continue
                sys.stdout.write('    angle_coeff @angle:'+angle_name+' '+
                                 angle2style[angle_name] +
                                 angle2params[angle_name] + 
                                 " # (ver=" + angle2ver[angle_name] +
                                 ", ref=" + angle2ref[angle_name] + ")\n")
                if angle_name in angle2class2_bb:
                    sys.stdout.write('    angle_coeff @angle:'+angle_name+' '+
                                     'bb ' + 
                                     angle2class2_bb[angle_name] +
                                     " # (ver=" + angle2ver[angle_name] +
                                     ", ref=" + angle2ref[angle_name] + ")\n")

                    assert(angle_name in angle2class2_ba)
                    sys.stdout.write('    angle_coeff @angle:'+angle_name+' '+
                                     'ba ' + 
                                     angle2class2_ba[angle_name] +
                                     " # (ver=" + angle2ver[angle_name] +
                                     ", ref=" + angle2ref[angle_name] + ")\n")
            sys.stdout.write('  }  # end of angle_coeff commands\n'
                             '\n\n')







        ################# Print 4-body Dihedral Interactions ##################

        dihedral_names_priority_high_to_low = sorted(dihedral2priority.items(),
                                                     key=itemgetter(1),
                                                     reverse=True)

        if len(dihedral_names_priority_high_to_low) > 0:
            sys.stdout.write('  # --------------- Dihedral Interactions: ---------------------\n')
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (4-body) "dihedral" interactions: --\n'
                             '  #  DihedralType AtmType1 AtmType2 AtmType3 AtmType3 [BondType1 Bnd2 Bnd3]\n')
            sys.stdout.write('\n\n'
                             '  write_once("Data Dihedrals By Type") {\n')
            for dihedral_name, in dihedral_names_priority_high_to_low:
                if not (dihedral2style[dihedral_name] in
                        dihedral_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(dihedral_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(dihedral_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @dihedral:' + dihedral_name + ' ' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade'
                                         + anames[0] +
                                         ',adc*,aie*,aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc'
                                         + anames[1] +
                                         ',aie*,aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc'
                                         + anames[2] +
                                         ',aie*,aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade'
                                         + anames[3] +
                                         ',adc*,aie*,aic*' +
                                         '\n')
                    else:
                        continue
                else:
                    sys.stdout.write('    @dihedral:' + dihedral_name + ' ' +
                                     ' @atom:*,b*,a*,d' + anames[0] + ',i* ' +
                                     ' @atom:*,b*,a*,d' + anames[1] + ',i* ' +
                                     '\n')

            sys.stdout.write('  }  # end of "Data Dihedrals By Type" section\n'
                             '\n')

            # Print the force-field parameters for these dihedral interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Dihedral Force Field Parameters: -------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write('    # '+dihedral_style2docs[dihedral_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # dihedral_coeff DihedralTypeName  DihedralStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for dihedral_name, in dihedral_names_priority_high_to_low:
                if not (dihedral2style[dihedral_name] in
                        dihedral_styles_selected):
                    continue
                # Did the user ask us to include "auto" interactions?
                if (IsAutoInteraction(dihedral_name) and
                    (not include_auto_equivalences)):
                    continue
                sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                 dihedral2style[dihedral_name] +
                                 dihedral2params[dihedral_name] +
                                 " # (ver=" + dihedral2ver[dihedral_name] +
                                 ", ref=" + dihedral2ref[dihedral_name] + ")\n")
                if dihedral_name in dihedral2class2_mbt:
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                     'mbt ' + 
                                     dihedral2class2_mbt[dihedral_name] +
                                     " # (ver=" + dihedral2ver[dihedral_name] +
                                     ", ref=" + dihedral2ref[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_ebt)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                     'ebt ' + 
                                     dihedral2class2_ebt[dihedral_name] +
                                     " # (ver=" + dihedral2ver[dihedral_name] +
                                     ", ref=" + dihedral2ref[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_at)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                     'at ' + 
                                     dihedral2class2_at[dihedral_name] +
                                     " # (ver=" + dihedral2ver[dihedral_name] +
                                     ", ref=" + dihedral2ref[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_aat)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                     'aat ' + 
                                     dihedral2class2_aat[dihedral_name] +
                                     " # (ver=" + dihedral2ver[dihedral_name] +
                                     ", ref=" + dihedral2ref[dihedral_name] + ")\n")
                    assert(dihedral_name in dihedral2class2_bb13)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dihedral_name+' '+
                                     'bb13 ' + 
                                     dihedral2class2_bb13[dihedral_name] +
                                     " # (ver=" + dihedral2ver[dihedral_name] +
                                     ", ref=" + dihedral2ref[dihedral_name] + ")\n")
            sys.stdout.write('  }  # end of dihedral_coeff commands\n'
                             '\n\n')





        ################# Print 4-body Improper Interactions ##################

        improper_names_priority_high_to_low = sorted(improper2priority.items(),
                                                     key=itemgetter(1),
                                                     reverse=True)

        if len(improper_names_priority_high_to_low) > 0:
            sys.stdout.write("  # --------------- Improper Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (4-body) "improper" interactions: --\n'
                             '  #  ImproperType AtmType1 AtmType2 AtmType3 AtmType3 [BondType1 Bnd2 Bnd3]\n')
            sys.stdout.write('\n'
                             '  write_once("Data Impropers By Type") {\n')
            for improper_name, in improper_names_priority_high_to_low:
                if not (improper2style[improper_name] in
                        improper_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(improper_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(improper_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @improper:' + improper_name +' '+
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc*,aie'
                                         + anames[0] + ',aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc*,aie*,aic'
                                         + anames[1] +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc*,aie'
                                         + anames[2] + ',aic*' +
                                         ' @atom:*,aq*,ab*,aae*,aac*,ade*,adc*,aie'
                                         + anames[3] + ',aic*' +
                                         '\n')
                    else:
                        continue
                else:
                    sys.stdout.write('    @improper:' + improper_name + ' ' +
                                     ' @atom:*,b*,a*,d*,i' + anames[0] + 
                                     ' @atom:*,b*,a*,d*,i' + anames[1] +
                                     ' @atom:*,b*,a*,d*,i' + anames[2] +
                                     ' @atom:*,b*,a*,d*,i' + anames[3] +
                                     '\n')

            sys.stdout.write('  }  # end of "Data Impropers By Type" section\n'
                             '\n')

            # Print the force-field parameters for these improper interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Improper Force Field Parameters: -------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write('    # '+improper_style2docs[improper_style]+'\n')
            sys.stdout.write('\n'
                             '# Syntax:  \n'
                             '  # improper_coeff ImproperTypeName  ImproperStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for improper_name, in improper_names_priority_high_to_low:
                if not (improper2style[improper_name] in
                        improper_styles_selected):
                    continue
                # Did the user ask us to include "auto" interactions?
                if (IsAutoInteraction(improper_name) and
                    (not include_auto_equivalences)):
                    continue
                sys.stdout.write('    improper_coeff @improper:'+improper_name+' '+
                                 improper2style[improper_name] +
                                 improper2params[improper_name] +
                                 " # (ver=" + improper2ver[improper_name] +
                                 ", ref=" + improper2ref[improper_name] + ")\n")
                if improper_name in improper2class2_aa:
                    sys.stdout.write('    improper_coeff @improper:'+improper_name+' '+
                                     'aa ' + 
                                     improper2class2_aa[improper_name] +
                                     " # (ver=" + improper2ver[improper_name] +
                                     ", ref=" + improper2ref[improper_name] + ")\n")
            sys.stdout.write('  }  # end of improper_coeff commands\n'
                             '\n\n')



        sys.stdout.write('\n\n\n\n'
                         '  # -------------------- Select LAMMPS style(s) ------------------\n'
                         '\n')

        
        sys.stdout.write('\n'
                         '  # LAMMPS supports many different kinds of bonded and non-bonded\n'
                         '  # interactions which can be selected at run time.  Eventually\n'
                         '  # we must inform LAMMPS which of them we will need.  We specify\n'
                         '  # this in the "In Init" section: \n\n')
        
        sys.stdout.write('  write_once("In Init") {\n')
        sys.stdout.write('    units real\n')
        sys.stdout.write('    atom_style full\n')

        if len(bond_styles) > 0:
            sys.stdout.write('    bond_style hybrid')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write(' ' + bond_style)
            sys.stdout.write('\n')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write('    # '+bond_style2docs[bond_style]+'\n')
            sys.stdout.write('\n')

        if len(angle_styles) > 0:
            sys.stdout.write('    angle_style hybrid')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write(' ' + angle_style)
            sys.stdout.write('\n')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write('    # '+angle_style2docs[angle_style]+'\n')
            sys.stdout.write('\n')

        if len(dihedral_styles) > 0:
            sys.stdout.write('    dihedral_style hybrid')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write(' ' + dihedral_style)
            sys.stdout.write('\n')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write('    # '+dihedral_style2docs[dihedral_style]+'\n')
            sys.stdout.write('\n')

        if len(improper_styles) > 0:
            sys.stdout.write('    improper_style hybrid')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write(' ' + improper_style)
            sys.stdout.write('\n')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write('    # '+improper_style2docs[improper_style]+'\n')
            sys.stdout.write('\n')

        if len(pair_styles) > 0:
            sys.stdout.write('    pair_style hybrid')
            for pair_style in pair_styles:
                if not (pair_style in pair_styles_selected):
                    continue
                sys.stdout.write(' ' + pair_style +
                                 ' ' + pair_style_args[pair_style])
            sys.stdout.write('\n')
            for pair_style in pair_styles:
                sys.stdout.write('    # '+pair_style2docs[pair_style]+'\n')
            sys.stdout.write('\n')

        sys.stdout.write('    pair_modify mix ' + pair_mixing_style + '\n')
        sys.stdout.write('    ' + special_bonds_command + '\n')
        sys.stdout.write('    ' + kspace_style + '\n')
        sys.stdout.write('  } #end of init parameters\n\n')
        sys.stdout.write('}  # ' + ffname + '\n\n')
        
        
        sys.stdout.write("#\n"
                         "# WARNING: The following 1-2, 1-3, and 1-4 weighting parameters were ASSUMED:\n")
        sys.stdout.write("#          " + special_bonds_command + "\n")
        sys.stdout.write(
            "#          (See http://lammps.sandia.gov/doc/special_bonds.html for details)\n")

        #sys.stderr.write(' done.\n')
        
        if filename_in != '':
            file_in.close()




    except InputError as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(1)



if __name__ == '__main__':
    main()
