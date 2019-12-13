#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2013

"""
ltemplify.py

The "ltemplify.py" script can be used to convert existing LAMMPS
input script and data files into a single .lt file
(which typically includes geometry, topology and force-field
 information for a single molecule in your system).

Example:

   ltemplify.py -name Mol file.in file.data > mol.lt

This creates a template for a new type of molecule (named "Mol"),
consisting of all the atoms in the lammps files you included,
and saves this data in a single ttree file ("mol.lt").
This file can be used with moltemplate to
define large systems containing this molecule.

"""

import sys
from collections import defaultdict

if sys.version < '2.6':
    raise InputError('Error: Using python ' + sys.version + '\n'
                     '       Alas, you must upgrade to a newer version of python (2.7 or later).')
elif sys.version > '3':
    from io import StringIO
else:
    from StringIO import StringIO

try:
    from .ttree_lex import *
    from .lttree_styles import *
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree_lex import *
    from lttree_styles import *

g_program_name = __file__.split('/')[-1]  # = 'ltemplify.py'
g_version_str = '0.65.0'
g_date_str = '2019-12-12'

def Intify(s):
    if s.isdigit():
        return int(s)
    elif s[0:2] == 'id':
        return int(s[2:])
    elif s[0:4] == 'type':
        return int(s[4:])
    else:
        return s


def IsNumber(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


def StringToInterval(sel_str, slice_delim='*'):
    # Split a string into 1-3 tokens using the slice_delim and convert to int.
    # What a mess. I should rewrite this function

    i_slice = sel_str.find(slice_delim)

    if i_slice == -1:
        a = sel_str
        b = sel_str
        c = ''
    else:
        a = sel_str[:i_slice]
        bc = sel_str[i_slice + len(slice_delim):]
        b = ''
        c = ''
        i_slice = bc.find(slice_delim)
        if i_slice == -1:
            b = bc
            c = ''
        else:
            b = bc[:i_slice]
            c = bc[i_slice + len(slice_delim):]

    if a == '':
        a = None
    elif a.isdigit():
        a = int(a)
    else:
        raise InputError('Error: invalid selection string \"' +
                         sel_str + '\"\n')

    if b == '':
        b = None
    elif b.isdigit():
        b = int(b)
    else:
        raise InputError('Error: invalid selection string \"' +
                         sel_str + '\"\n')

    if c == '':
        c = None
    elif c.isdigit():
        c = int(c)
    else:
        raise InputError('Error: invalid selection string \"' +
                         sel_str + '\"\n')

    if c == None:
        return (a, b)
    else:
        return (a, b, c)


# Selections are simply lists of 2-tuples (pairs)

def LammpsSelectToIntervals(sel_str, slice_delim='*', or_delim=', '):
    """
    This function converts a string such as "1*4 6 9*12 50*70*10" into
    a list of tuples, for example: [(1,4), (6,6), (9,12), (50,50), (60,60), (70,70)]
    In general, the of intervals has the form:
       [(a1,b1), (a2,b2), (a3,b3), ... ]

    An atom is considered to belong to this selection
    if it happens to lie within the closed interval [a,b]
    for any pair of a,b values in the list of intervals.
    If for a given pair a,b, either a or b is "None", then that a or b
    value is not used to disqualify membership in the interval.
    (Similar to -infinity or +infinity.  In other words if a is set to None,
     then to belong to the interval it is enough to be less than b.)

    """
    selection_list = []
    # tokens = sel_str.split(or_delim) <-- Not what we want when
    # len(or_delim)>1
    tokens = LineLex.TextBlock2Lines(sel_str, or_delim, keep_delim=False)
    for token in tokens:
        token = token.strip()
        interval = StringToInterval(token, slice_delim)

        if len(interval) == 2:
            # Normally, "interval" should be a tuple containing 2 entries
            selection_list.append(interval)
        else:
            assert(len(interval) == 3)
            # Handle 1000:2000:10 notation
            # (corresponding to 1000, 1010, 1020, 1030, ..., 1990, 2000)
            a = interval[0]
            b = interval[1]
            incr = interval[2]
            i = a
            while i <= b:
                selection_list.append((i, i))
                i += incr

    return selection_list


def IntervalListToMinMax(interval_list):
    min_a = None
    max_b = None
    for (a, b) in interval_list:
        if ((not (type(a) is int)) or (not (type(b) is int))):
            return None, None  # only integer min/max makes sense. otherwise skip

        if (min_a == None) or (a < min_a):
            min_a = a
        if (max_b == None) or (b > max_b):
            max_b = b
    return min_a, max_b


def MergeIntervals(interval_list):
    """
    A crude simple function that merges consecutive intervals in the list
    whenever they overlap.  (This function does not bother to compare
    non-consecutive entries in the interval_list.)

    """
    i = 1
    while i < len(interval_list):
        if ((interval_list[i - 1][1] == None) or
                (interval_list[i - 1][1] + 1 >= interval_list[i][0])):
            interval_list[i - 1] = (interval_list[i - 1]
                                    [0], interval_list[i][1])
            del interval_list[i]
        else:
            i += 1


def BelongsToSel(i, sel):
    if (i == None) or (sel == None) or (len(sel) == 0):
        # If the user has not specified a selection for this category,
        # then by default all objects are accepted
        return True

    elif (type(i) is str):
        if i.isdigit():
            i = int(i)
        else:
            return True

    belongs = False
    for interval in sel:
        assert(len(interval) == 2)
        if interval[0]:
            if i >= interval[0]:
                if (interval[1] == None) or (i <= interval[1]):
                    belongs = True
                    break
        elif interval[1]:
            if i <= interval[1]:
                belongs = True
                break
        else:
            # In that case, the user entered something like "*"
            # which covers all possible numbers
            belongs = True
            break

    return belongs



def LookupVarName(i, int2name=None, default_name = ''):
    if int2name and (i in int2name):
        var_name = int2name[i]
    else:
        var_name = default_name + str(i)
    return var_name



def Stringify(i, int2name, prefix1, prefix2, default_name = ''):
    var_name = LookupVarName(i, int2name, default_name)
    if var_name.find(' ') != -1:
        formatted_var_name = prefix1 + '{' + prefix2 + ':' + var_name + '}'
    else:
        formatted_var_name = prefix1 + prefix2 + ':' + var_name
    return formatted_var_name



class Ltemplify(object):

    def __init__(self, argv):
        """
        The constructor requires a list of command line arguments to 
        figure out the format of the output file we will generate.
        This meaning of these arguments is explained in "doc_ltemplify.md":
       https://github.com/jewettaij/moltemplate/blob/master/doc/doc_ltemplify.md
        Note: You can either specify the input scripts and data files in the
        argument list, OR specify them later by passing them as arguments
        to the Convert() member function.  The second approach is preferred.
        """

        self.atom_column_names = []
        self.non_empty_output = False
        self.no_warnings = True
        self.indent = 2
        self.cindent = 0
        self.atomid_selection = []
        self.atomtype_selection = []
        self.molid_selection = []
        self.mol_name = ''

        self.min_sel_atomid = None
        self.min_sel_atomtype = None
        self.min_sel_bondid = None
        self.min_sel_bondtype = None
        self.min_sel_angleid = None
        self.min_sel_angletype = None
        self.min_sel_dihedralid = None
        self.min_sel_dihedraltype = None
        self.min_sel_improperid = None
        self.min_sel_impropertype = None

        self.max_sel_atomid = None
        self.max_sel_atomtype = None
        self.max_sel_bondid = None
        self.max_sel_bondtype = None
        self.max_sel_angleid = None
        self.max_sel_angletype = None
        self.max_sel_dihedralid = None
        self.max_sel_dihedraltype = None
        self.max_sel_improperid = None
        self.max_sel_impropertype = None

        self.needed_atomids = set([])
        self.needed_atomtypes = set([])
        self.needed_molids = set([])
        self.needed_bondids = set([])
        self.needed_bondtypes = set([])
        self.needed_angleids = set([])
        self.needed_angletypes = set([])
        self.needed_dihedralids = set([])
        self.needed_dihedraltypes = set([])
        self.needed_improperids = set([])
        self.needed_impropertypes = set([])
        self.needed_cmapids = set([])
        self.needed_cmaptypes = set([])

        self.min_needed_atomtype = None
        self.max_needed_atomtype = None
        self.min_needed_bondtype = None
        self.max_needed_bondtype = None
        self.min_needed_angletype = None
        self.max_needed_angletype = None
        self.min_needed_dihedraltype = None
        self.max_needed_dihedraltype = None
        self.min_needed_impropertype = None
        self.max_needed_impropertype = None
        self.min_needed_cmaptype = None
        self.max_needed_cmaptype = None

        self.min_needed_atomid = None
        self.max_needed_atomid = None
        self.min_needed_molid = None
        self.max_needed_molid = None
        self.min_needed_bondid = None
        self.max_needed_bondid = None
        self.min_needed_angleid = None
        self.max_needed_angleid = None
        self.min_needed_dihedralid = None
        self.max_needed_dihedralid = None
        self.min_needed_improperid = None
        self.max_needed_improperid = None
        self.min_needed_cmapid = None
        self.max_needed_cmapid = None

        # To process the selections, we need to know the atom style:
        self.atom_style_undefined = True

        self.i_atomid = None
        self.i_atomtype = None
        self.i_molid = None
        self.i_x = None
        self.i_y = None
        self.i_z = None

        self.l_in_init = []
        self.l_in_settings = []
        self.l_in_masses = []
        self.l_in_pair_coeffs = []
        self.l_in_bond_coeffs = []
        self.l_in_angle_coeffs = []
        self.l_in_dihedral_coeffs = []
        self.l_in_improper_coeffs = []
        self.l_in_group = []
        self.l_in_set = []
        self.l_in_set_static = []
        self.l_in_fix_shake = []
        self.l_in_fix_rigid = []
        self.l_in_fix_poems = []
        self.l_in_fix_qeq = []
        self.l_in_fix_qmmm = []
        self.l_data_masses = []
        self.l_data_bond_coeffs = []
        self.l_data_angle_coeffs = []
        self.l_data_dihedral_coeffs = []
        self.l_data_improper_coeffs = []
        self.l_data_cmap_coeffs = []
        self.l_data_pair_coeffs = []
        self.l_data_pairij_coeffs = []
        self.l_data_atoms = []
        self.l_data_velocities = []
        self.l_data_bonds = []
        self.l_data_angles = []
        self.l_data_dihedrals = []
        self.l_data_impropers = []
        self.l_data_cmap = []

        # class2 force fields
        # l_in_bondbond_coeffs = []   <--not needed, included in l_in_angle_coeff
        # l_in_bondangle_coeffs = []   <--not needed, included in l_in_angle_coeff
        # l_in_middlebondtorsion_coeffs = [] not needed, included in l_in_dihedral_coeff
        # l_in_endbondtorsion_coeffs = [] <--not needed, included in l_in_dihedral_coeff
        # l_in_angletorsion_coeffs = [] <--not needed, included in l_in_dihedral_coeff
        # l_in_angleangletorsion_coeffs = [] not needed, included in l_in_dihedral_coeff
        # l_in_bondbond13_coeffs = []  <--not needed, included in l_in_dihedral_coeff
        # l_in_angleangle_coeffs = []  <--not needed, included in
        # l_in_improper_coeff
        self.l_data_bondbond_coeffs = []
        self.l_data_bondangle_coeffs = []
        self.l_data_middlebondtorsion_coeffs = []
        self.l_data_endbondtorsion_coeffs = []
        self.l_data_angletorsion_coeffs = []
        self.l_data_angleangletorsion_coeffs = []
        self.l_data_bondbond13_coeffs = []
        self.l_data_angleangle_coeffs = []

        # non-point-like particles:
        self.l_data_ellipsoids = []
        self.l_data_lines = []
        self.l_data_triangles = []

        # automatic generation of bonded interactions by type:
        self.l_data_angles_by_type = []
        self.l_data_dihedrals_by_type = []
        self.l_data_impropers_by_type = []

        self.atoms_already_read = False
        self.some_pair_coeffs_read = False
        self.complained_atom_style_mismatch = False
        self.infer_types_from_comments = True
        self.ignore_coeffs = False
        self.ignore_angles_dihedrals_impropers = False
        self.ignore_bond_types = False
        self.ignore_masses = False
        self.forbid_type_name_duplicates = False
        self.prepend_atom_type = ''
        self.remove_coeffs_from_data_file = True


        # self.input_data_file      Name of the data file we will read.
        # self.input_script_files   Name of the LAMMPS input scripts to be read.
        self.input_data_file = None
        self.input_script_files = []
        # Normally these self.input_file... data members are not necessary
        # when this function is invoked from within python.
        # Python users will specify the list of files they want to read
        # by passing arguments to the "Convert()" member function.
        # These members are only needed when this function is used to parse
        # all of the arguments from the command line.  Those arguments will
        # include a list of the names of input files that we need to read.
        # The next two data members are then used to store that information.


        # Parse the argument list.
        # Arguments are explained in the ltemplify.py documentation located at:
        # https://github.com/jewettaij/moltemplate/blob/master/doc/utils/doc_ltemplify.md

        i = 0

        while i < len(argv):

            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

            if argv[i] == '-columns':
                if i + 1 >= len(argv):
                    raise InputError('Error: the \"' + argv[i] + '\" argument should be followed by a quoted\n'
                                     '       string which contains a space-delimited list of the names of\n'
                                     '       of columns in the \"Atoms\" section of the LAMMPS data file.\n'
                                     '       If the list contains the symbols:\n'
                                     '    \"atom-ID\" or \"atomid\", they are interpreted\n'
                                     '       as unique atom ID numbers, and columns named\n'
                                     '    \"atom-type\" or \"atomtype\" are interpreted\n'
                                     '       as atom types.  Finally, columns named\n'
                                     '    \"molecule-ID\", \"molecule\", or \"mol-ID\", or \"mol\"\n'
                                     '       are interpreted as unique molecule id numbers.\n'
                                     'Example:\n'
                                     '    ' +
                                     argv[
                                         i] + ' \'atom-ID atom-type q polarizability molecule-ID x y z\'\n'
                                     '       defines a custom atom_style containing the properties\n'
                                     '            atom-ID atom-type q polarizability molecule-ID x y z\n'
                                     '    Make sure you enclose the entire list in quotes.\n')
                self.atom_column_names = argv[i+1].strip('\"\'').strip().split()
                del argv[i:i + 2]

            elif ((argv[i] == '-name') or
                  (argv[i] == '-molname') or
                  (argv[i] == '-molecule-name') or
                  (argv[i] == '-molecule_name')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a a molecule type name.\n')
                self.cindent = 2
                self.indent += self.cindent
                self.mol_name = argv[i + 1]
                del argv[i:i + 2]

            elif ((argv[i].lower() == '-atomstyle') or
                  (argv[i].lower() == '-atom_style') or
                  (argv[i].lower() == '-atom-style')):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a an atom_style name.\n'
                                     '       (or single quoted string which includes a space-separated\n'
                                     '       list of column names).\n')
                self.atom_style_undefined = False
                self.atom_column_names = AtomStyle2ColNames(argv[i + 1])
                if (argv[i + 1].strip().split()[0] in g_style_map):
                    self.l_in_init.append((' ' * self.indent) +
                                     'atom_style ' + argv[i + 1] + '\n')
                sys.stderr.write('\n    \"Atoms\" column format:\n')
                sys.stderr.write('    ' + (' '.join(self.atom_column_names)) + '\n')
                self.i_atomid, self.i_atomtype, self.i_molid = ColNames2AidAtypeMolid(
                    self.atom_column_names)
                # Which columns contain the coordinates?
                ii_coords = ColNames2Coords(self.atom_column_names)
                assert(len(ii_coords) == 1)
                self.i_x = ii_coords[0][0]
                self.i_y = ii_coords[0][1]
                self.i_z = ii_coords[0][2]

                if self.i_molid:
                    sys.stderr.write('      (i_atomid=' + str(self.i_atomid + 1) + ', i_atomtype=' + str(
                        self.i_atomtype + 1) + ', i_molid=' + str(self.i_molid + 1) + ')\n\n')
                else:
                    sys.stderr.write('      (i_atomid=' + str(self.i_atomid + 1) +
                                     ', i_atomtype=' + str(self.i_atomtype + 1) + ')\n')
                del argv[i:i + 2]

            elif ((argv[i].lower() == '-id') or
                  (argv[i].lower() == '-atomid') or
                  #(argv[i].lower() == '-atomids') or
                  (argv[i].lower() == '-atom-id')
                   #(argv[i].lower() == '-atom-ids') or
                  #(argv[i].lower() == '-$atom') or
                  #(argv[i].lower() == '-$atoms')
                  ):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a list of integers\n'
                                     '       (or strings).  These identify the group of atoms you want to\n'
                                     '       to include in the template you are creating.\n')
                self.atomid_selection += LammpsSelectToIntervals(argv[i + 1])
                self.min_sel_atomid, self.max_sel_atomid = IntervalListToMinMax(
                    self.atomid_selection)
                del argv[i:i + 2]
            elif ((argv[i].lower() == '-datacoeffs') or
                  (argv[i].lower() == '-datacoeff') or
                  (argv[i].lower() == '-Coeff') or
                  (argv[i].lower() == '-Coeffs')):
                self.remove_coeffs_from_data_file = False
                del argv[i:i + 1]
            elif ((argv[i].lower() == '-type') or
                  #(argv[i].lower() == '-t') or
                  (argv[i].lower() == '-atomtype') or
                  (argv[i].lower() == '-atom-type')
                  #(argv[i].lower() == '-atomtypes') or
                  #(argv[i].lower() == '-atom-types') or
                  #(argv[i].lower() == '-@atom') or
                  #(argv[i].lower() == '-@atoms') or
                  #(argv[i].lower() == '-@atomtype') or
                  #(argv[i].lower() == '-@atomtypes')
                  ):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a list of integers.\n'
                                     '       (or strings).  These identify the group of atom types you want to\n'
                                     '       to include in the template you are creating.\n')
                self.atomtype_selection += LammpsSelectToIntervals(argv[i + 1])
                self.min_sel_atomtype, self.max_sel_atomtype = IntervalListToMinMax(
                    self.atomtype_selection)
                del argv[i:i + 2]
            elif ((argv[i].lower() == '-mol') or
                  #(argv[i].lower() == '-m') or
                  (argv[i].lower() == '-molid') or
                  #(argv[i].lower() == '-molids') or
                  (argv[i].lower() == '-mol-id') or
                  #(argv[i].lower() == '-mol-ids') or
                  #(argv[i].lower() == '-molecule') or
                  (argv[i].lower() == '-moleculeid') or
                  (argv[i].lower() == '-molecule-id')
                  #(argv[i].lower() == '-molecules') or
                  #(argv[i].lower() == '-molecule-ids') or
                  #(argv[i].lower() == '-$mol') or
                  #(argv[i].lower() == '-$molecule')
                  ):
                if i + 1 >= len(argv):
                    sys.stderr.write('Error: ' + argv[i] + ' flag should be followed by a list of integers.\n'
                                     '       (or strings).  These identify the group of molecules you want to\n'
                                     '       include in the template you are creating.\n')
                self.molid_selection += LammpsSelectToIntervals(argv[i + 1])
                del argv[i:i + 2]

            elif (argv[i] == '-ignore-comments'):
                self.infer_types_from_comments = False
                del argv[i:i + 1]

            elif (argv[i] == '-infer-comments'):
                self.infer_types_from_comments = True
                del argv[i:i + 1]

            elif (argv[i] == '-ignore-coeffs'):
                self.ignore_coeffs = True
                del argv[i:i + 1]

            elif (argv[i] == '-ignore-angles'):
                self.ignore_angles_dihedrals_impropers = True
                del argv[i:i + 1]

            elif (argv[i] == '-ignore-bond-types'):
                self.ignore_bond_types = True
                del argv[i:i + 1]

            elif (argv[i] == '-ignore-masses'):
                self.ignore_masses = True
                del argv[i:i + 1]

            elif (argv[i] == '-ignore-duplicates'):
                self.forbid_type_name_duplicates = False
                del argv[i:i + 1]

            elif (argv[i] == '-forbid-duplicates'):
                self.forbid_type_name_duplicates = True
                del argv[i:i + 1]

            elif (argv[i] == '-prepend-atom-type'):
                self.prepend_atom_type = argv[i+1]
                if ((len(self.prepend_atom_type) > 0) and
                    (self.prepend_atom_type[-1] != '/')):
                    self.prepend_atom_type += '/'
                del argv[i:i + 2]

            else:
                i += 1

        # We might need to parse the simulation boundary-box.
        # If so, use these variables.  (None means uninitialized.)
        self.boundary_xlo = None
        self.boundary_xhi = None
        self.boundary_ylo = None
        self.boundary_yhi = None
        self.boundary_zlo = None
        self.boundary_zhi = None
        self.boundary_xy = None
        self.boundary_yz = None
        self.boundary_xz = None

        # atom type names
        #self.atomtypes_name2int = {}
        self.atomtypes_name2list = defaultdict(list)
        self.atomtypes_name2int = {}
        self.atomtypes_int2name = {}
        # self.atomids_name2int = {}  not needed
        self.atomids_int2name = {}
        self.atomids_by_type = defaultdict(list)

        self.atomid2type = {}
        self.atomid2mol = {}

        # molecule id names
        self.molids_int2name = {} #(empty for now. perhaps I'll fill this later)

        # bond type names
        self.bondtypes_int2name = {}

        # angle type names
        self.angletypes_int2name = {}

        # dihedral type names
        self.dihtypes_int2name = {}

        # improper type names
        self.imptypes_int2name = {}

        #---------------------------------------------------------
        #-- The remaining arguments are files that the user wants
        #-- us to read and convert.  It is typical to have
        #-- multiple input files, because LAMMPS users often
        #-- store their force field parameters in either the LAMMPS
        #-- data files and input script files, or both.
        #-- We want to search all of the LAMMPS input files in
        #-- order to make sure we extracted all the force field
        #-- parameters (coeff commands).
        #---------------------------------------------------------

        if len(argv) > 0:
            self.input_data_file = argv[-1]  #the last argument is the data file
            self.input_script_files = argv[:-1] # optional input script files






    def Convert(self,
                out_file,
                input_data_file=None,
                input_script_files=None):

        """
        Converts a data file (and, optionally, one or more input script files)
        into a new file ("out_file") which is in MOLTEMPLATE (.LT) format.
        The arguments can be either filenames or StringIO objects.
        The "input_script_file" argument can be a single string or StringIO
        object, or a list of such objects.
        (See doc_ltemplify.md for details.)
        """
        # We need to know the list of files we will read.

        # If the caller did not specify parameters, then copy them
        # from self.input_data_file and self.input_script_files.
        # If the caller did specify these files, this will override the list of
        # files stored in in self.input_data_file and self.input_script_files.
        if input_data_file == None:
            input_data_file = self.input_data_file
        assert(input_data_file != None)
        if input_script_files == None:
            input_script_files = self.input_script_files
        elif (not hasattr(input_script_files, '__getitem__')):
            input_script_files = [input_script_files]

        # (Note: The "data" file is assumed to be the last entry in the list.)
        input_files = input_script_files + [input_data_file]
        assert(len(input_files) > 0)
        input_files_txt = ['' for f in input_files]
        input_files_sio = [None for f in input_files]

        #### READ ALL OF THE INPUT FILES AND SAVE THEIR TEXT FOR LATER
        for i_f in range(0, len(input_files)):
            fname = input_files[i_f]
            if isinstance(fname, str):
                try:
                    input_file = open(fname, 'r')
                except IOError:
                    raise InputError('Error: unrecognized argument (\"' + fname + '\"),\n'
                                     '       OR unable to open file:\n'
                                     '\n'
                                     '       \"' + fname + '\"\n'
                                     '       for reading.\n'
                                     '\n'
                                     '       (If you were not trying to open a file with this name,\n'
                                     '        then there is a problem in your argument list.)\n')

                sys.stderr.write('reading file \"' + fname + '\"\n')
            else:
                input_file = fname #then "fname" is a file stream, not a string

            input_files_txt[i_f] = input_file.read()
        #### FINISHED READING THE FILES

        # Now create StringIO versions of these files
        for i_f in range(0, len(input_files)):
            input_files_sio[i_f] = StringIO(input_files_txt[i_f])

        # Split this list of files into a LAMMPS data files and input scripts:
        # (Note: The "data" file is assumed to be the last entry in the list.)
        assert(len(input_files_sio) > 0)
        input_data_file_sio = input_files_sio[-1]
        input_script_files_sio = input_files_sio[:-1]


        # PASS 1
        # Determine the atom_style, as well as the atom type names
        self.Pass1(input_data_file_sio, input_script_files_sio)


        # PASS 2
        # Parse all other sections of the LAMMPS files,
        # including Atoms, Bonds, Angles, Dihedrals, Impropers and Masses

        # Again, create StringIO versions of these files for reading again.
        for i_f in range(0, len(input_files)):
            input_files_sio[i_f] = StringIO(input_files_txt[i_f])
        input_data_file_sio = input_files_sio[-1]
        input_script_files_sio = input_files_sio[:-1]

        self.Pass2(input_data_file_sio, input_script_files_sio)



        # POST PROCESSING:
        # Deal with CUSTOM atom type and bonded interaction type names.
        #
        # Infer the atomtype names, (and bondtype names, angletype names,
        # dihedraltype names, and impropertype names)
        # ...from comments located in the data file
        self.PostProcess1()

        # Determine atomid names from atomtype names, and update the references
        # to these atomids the "Atoms", "Velocities", "Ellipsoids",... sections.
        self.PostProcess2()

        # Discard unnecessary information.
        # Delete the bonded and nonbonded interactions involving atoms
        # that we don't care about (atoms that were not selected by the user).
        # This is also a good time to update the references to the atomids in
        # bonded interactions (to refer to custom atomid names, if applicable).
        self.PostProcess3()

        # Print all the information contained in the list variables
        # (such as "l_in_..." and "l_data_...") to a file (out_file).
        self.Write(out_file)



    def Pass1(self, input_data_file, input_script_files):
        "PASS1: determine the atom_style, as well as the atom type names."

        sys.stderr.write(Ltemplify.Pass1.__doc__+'\n')

        atom_style_str = ''

        input_files = input_script_files + [input_data_file]

        for i_arg in range(0, len(input_files)):
            fname = input_files[i_arg]
            if isinstance(fname, str):
                try:
                    lammps_file = open(fname, 'r')
                except IOError:
                    raise InputError('Error: unrecognized argument (\"' + fname + '\"),\n'
                                     '       OR unable to open file:\n'
                                     '\n'
                                     '       \"' + fname + '\"\n'
                                     '       for reading.\n'
                                     '\n'
                                     '       (If you were not trying to open a file with this name,\n'
                                     '        then there is a problem in your argument list.)\n')

                sys.stderr.write('reading file \"' + fname + '\"\n')
            else:
                lammps_file = fname #then "fname" is a file stream, not a string


            lex = LineLex(lammps_file, fname)
            lex.source_triggers = set(['include', 'import'])
            lex.commenters = '' # do not ignore comments
            # set up lex to accept most characters in file names:
            lex.wordterminators = '(){}' + lex.whitespace
            # set up lex to understand the "include" statement:
            lex.source = 'include'
            lex.escape = '\\'

            while lex:

                line_orig = lex.ReadLine()
                line = line_orig
                comment_text = ''
                ic = line_orig.find('#')
                if ic != -1:
                    comment_text = line[ic + 1:].strip()
                    line = line_orig[:ic]

                line = line.strip()

                #sys.stderr.write('  processing \"'+line.strip()+'\", (\"'+lex.infile+'\":'+str(lex.lineno)+')\n')

                tokens = line.split()
                if (len(tokens) > 0):
                    if ((tokens[0] == 'atom_style') and
                        self.atom_style_undefined):

                        if self.atom_style_undefined:
                            self.atom_style_undefined = False
                            atom_style_str = comment_text.strip()

                        sys.stderr.write(
                            '  Atom Style found. Processing: \"' + line + '\"\n')
                        if self.atoms_already_read:
                            raise InputError('Error: The file containing the \"atom_style\" command must\n'
                                             '       come before the data file in the argument list.\n'
                                             '       (The templify program needs to know the atom style before reading\n'
                                             '       the data file.  Either change the order of arguments so that the\n'
                                             '       LAMMPS input script file is processed before the data file, or use\n'
                                             '       the \"-atom_style\" command line argument to specify the atom_style.)\n')

                        atom_style_str = ' '.join(tokens[1:])  # skip over the 'atom_style '

                    elif (line == 'Atoms'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        # Check to see if the atom_style was specified
                        # in a comment at the end of this line.

                        if self.atom_style_undefined:
                            self.atom_style_undefined = False
                            atom_style_str = comment_text.strip()


                    elif (line == 'Masses'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            # Read the next line of text but don't skip comments
                            #comment_char_backup = lex.commenters
                            lex.commenters = ''
                            line_orig = lex.ReadLine()
                            #lex.commenters = comment_char_backup

                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.rstrip()

                            if line.strip() in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.strip().split()
                            if len(tokens) > 0:
                                atomtype = Intify(tokens[0])
                                atomtype_name = 'type' + str(atomtype)

                                if comment_text != '':

                                    # Assume the entire comment is the atom type
                                    atomtype_name = comment_text.strip()

                                    # COMMENTING OUT:
                                    ## Assume the first word after the comment
                                    ## character is the atom type name
                                    #comment_tokens = comment_text.split()
                                    ##comment_tokens = SplitQuotedString(comment_text,
                                    ##                                   comment_char='')
                                    #atomtype_name = (self.prepend_atom_type +
                                    #                 comment_tokens[0])

                                if BelongsToSel(atomtype, self.atomtype_selection):

                                    #Infer atom type names from comment strings?
                                    if self.infer_types_from_comments:
                                        self.atomtypes_name2list[atomtype_name].append(atomtype)
                                    else:
                                        self.atomtypes_int2name[atomtype] = ('type' +
                                                                             str(atomtype))

                        if self.infer_types_from_comments:
                            for atypename, ilist in self.atomtypes_name2list.items():
                                assert(len(ilist) > 0)
                                if len(ilist) == 1:
                                    self.atomtypes_int2name[ilist[0]] = atypename
                                    self.atomtypes_name2int[atypename] = ilist[0]
                            for atypename, ilist in self.atomtypes_name2list.items():
                                if len(ilist) > 1:
                                    if self.forbid_type_name_duplicates:
                                        raise InputError('Error: duplicate atom type names in Masses section: \"' + atypename + '\"\n'
                                                         '       (By default ' + g_program_name +
                                                         ' attempts to infer atom type names from\n'
                                                         '       comments which appear in the \"Masses\" section of your data file.)\n'
                                                         '       This error message occurs if two different atom types\n'
                                                         '       have the same name (in the comments section)\n'
                                                         '       You can avoid this error by adding the \"-ignore-duplicates\" argument.\n')
                                        
                                    # If there are multiple different atom types
                                    # corresponding to the same name (eg "C")
                                    # then add a numeric suffix to the name
                                    # (eg "C1", "C2", ...).
                                    ilist = sorted(ilist)
                                    isuffix = 1
                                    I = 0
                                    while I < len(ilist):
                                        name_attempt = atypename+str(isuffix)
                                        if ((name_attempt in self.atomtypes_name2int)
                                            or
                                            (name_attempt in self.atomtypes_name2list)):
                                            # Be careful not to create a name
                                            # that has already been used.
                                            # Keep incrementing "isuffix" until
                                            # "name_attempt" has not been used
                                            isuffix += 1
                                        else:
                                            self.atomtypes_name2int[name_attempt] = ilist[I]
                                            self.atomtypes_int2name[ilist[I]] = name_attempt
                                            I += 1

            # We are finished reading that file.  close it.
            lammps_file.close()

            # (As a C++ programmer, this is why I don't like the fact
            #  that python  uses indentation exclusively to indicate scope.
            #  To be fair, I should use more function calls...)

        # Loop over input files...

        # Now we are through reading all of the LAMMPS input files.

        if atom_style_str == '':
            # The default atom_style is "full"
            atom_style_str = 'full'

        sys.stderr.write('\n    atom_style = \"'+atom_style_str+'\"\n')

        self.atom_column_names = AtomStyle2ColNames(atom_style_str)
        self.i_atomid, self.i_atomtype, self.i_molid = ColNames2AidAtypeMolid(
            self.atom_column_names)
        # Which columns contain the coordinates?
        ii_coords = ColNames2Coords(self.atom_column_names)
        assert(len(ii_coords) == 1)
        self.i_x = ii_coords[0][0]
        self.i_y = ii_coords[0][1]
        self.i_z = ii_coords[0][2]

        sys.stderr.write('\n    \"Atoms\" column format:\n')
        sys.stderr.write('    ' + (' '.join(self.atom_column_names)) + '\n')
        if self.i_molid:
            sys.stderr.write('        (i_atomid=' +
                             str(self.i_atomid + 1) +
                             ', i_atomtype=' + str(
                                 self.i_atomtype + 1) +
                             ', i_molid=' +
                             str(self.i_molid + 1) + ')\n\n')
        else:
            sys.stderr.write('        (i_atomid=' +
                             str(self.i_atomid + 1) +
                             ', i_atomtype=' +
                             str(self.i_atomtype + 1) + ')\n\n')
        self.l_in_init.append((' ' * self.indent) +
                         'atom_style ' + atom_style_str)







    def Pass2(self, input_data_file, input_script_files):
        "PASS2: Parse Atoms, Bonds, Angles, Dihedrals, Impropers and Masses."

        sys.stderr.write(Ltemplify.Pass2.__doc__+'\n')

        input_files = input_script_files + [input_data_file]

        for i_arg in range(0, len(input_files)):
            fname = input_files[i_arg]
            if isinstance(fname, str):
                try:
                    lammps_file = open(fname, 'r')
                except IOError:
                    raise InputError('Error: unrecognized argument (\"' + fname + '\"),\n'
                                     '       OR unable to open file:\n'
                                     '\n'
                                     '       \"' + fname + '\"\n'
                                     '       for reading.\n'
                                     '\n'
                                     '       (If you were not trying to open a file with this name,\n'
                                     '        then there is a problem in your argument list.)\n')

                sys.stderr.write('reading file \"' + fname + '\"\n')
            else:
                lammps_file = fname #then "fname" is a file stream, not a string

            lex = LineLex(lammps_file, fname)
            lex.source_triggers = set(['include', 'import'])
            lex.commenters = '' # do not ignore comments
            # set up lex to accept most characters in file names:
            lex.wordterminators = '(){}' + lex.whitespace
            # set up lex to understand the "include" statement:
            lex.source = 'include'
            lex.escape = '\\'

            while lex:
                line_orig = lex.ReadLine()
                line = line_orig
                comment_text = ''
                ic = line.find('#')
                if ic != -1:
                    comment_text = line[ic + 1:].strip()
                    line = line[:ic]
                line = line.strip()

                #sys.stderr.write('  processing \"'+line+'\", (\"'+lex.infile+'\":'+str(lex.lineno)+')\n')

                tokens = line.split()
                if (len(tokens) > 0):

                    if (tokens[0] in set(['units',
                                          'angle_style',
                                          'bond_style',
                                          'dihedral_style',
                                          'improper_style',
                                          'min_style',
                                          'pair_style',
                                          'pair_modify',
                                          'special_bonds',
                                          'kspace_style',
                                          'kspace_modify'])):
                        self.l_in_init.append((' ' * self.indent) + line.lstrip())

                    # if (line == 'LAMMPS Description'):
                    #    sys.stderr.write('  reading \"'+line+'\"\n')
                    #    # skip over this section
                    #    while lex:
                    #        line = lex.ReadLine()
                    #        if line in data_file_header_names:
                    #            lex.push_raw_text(line+'\n') # <- Save line for later
                    #            break

                    elif (line == 'Atoms'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        self.atoms_already_read = True

                        # Before attempting to read atomic coordinates, first find
                        # the lattice vectors of the simulation's boundary box:
                        #    Why do we care about the Simulation Boundary?
                        # Some LAMMPS data files store atomic coordinates in a
                        # complex format with 6 numbers, 3 floats, and 3 integers.
                        # The 3 floats are x,y,z coordinates. Any additional numbers
                        # following these are integers which tell LAMMPS which cell
                        # the particle belongs to, (in case it has wandered out of
                        # the original periodic boundary box). In order to find
                        # the true location of the particle, we need to offset that
                        # particle's position with the unit-cell lattice vectors:
                        # avec, bvec, cvec  (or multiples thereof)
                        # avec, bvec, cvec are the axis of the parallelepiped which
                        # define the simulation's boundary. They are described here:
                        # http://lammps.sandia.gov/doc/Section_howto.html#howto-12
                        if ((self.boundary_xlo == None) or (self.boundary_xhi == None) or
                            (self.boundary_ylo == None) or (self.boundary_yhi == None) or
                            (self.boundary_zlo == None) or (self.boundary_zhi == None)):

                            raise InputError('Error: Either DATA file lacks a boundary-box header, or it is in the wrong\n'
                                             '       place.  At the beginning of the file, you need to specify the box size:\n'
                                             '       xlo xhi ylo yhi zlo zhi   (and xy xz yz if triclinic)\n'
                                             '       These numbers should appear BEFORE the other sections in the data file\n'
                                             '       (such as the \"Atoms\", \"Masses\", \"Bonds\", \"Pair Coeffs\" sections)\n'
                                             '\n'
                                             '       Use this format (example):\n'
                                             '       -100.0 100.0 xhi xlo\n'
                                             '        0.0  200.0  yhi ylo\n'
                                             '       -25.0 50.0   zhi zlo\n'
                                             '\n'
                                             'For details, see http://lammps.sandia.gov/doc/read_data.html\n'
                                             '\n'
                                             '       (NOTE: If the atom coordinates are NOT followed by integers, then\n'
                                             '       these numbers are all ignored, however you must still specify\n'
                                             '       xlo, xhi, ylo, yhi, zlo, zhi.  You can set them all to 0.0.)\n')

                        if not (self.boundary_xy and self.boundary_yz and self.boundary_xz):
                            # Then use a simple rectangular boundary box:
                            avec = (self.boundary_xhi - self.boundary_xlo, 0.0, 0.0)
                            bvec = (0.0, self.boundary_yhi - self.boundary_ylo, 0.0)
                            cvec = (0.0, 0.0, self.boundary_zhi - self.boundary_zlo)
                        else:
                            # Triclinic geometry in LAMMPS is explained here:
                            # http://lammps.sandia.gov/doc/Section_howto.html#howto-12
                            # http://lammps.sandia.gov/doc/read_data.html
                            avec = (self.boundary_xhi - self.boundary_xlo, 0.0, 0.0)
                            bvec = (self.boundary_xy, self.boundary_yhi - self.boundary_ylo, 0.0)
                            cvec = (self.boundary_xz, self.boundary_yz,
                                    self.boundary_zhi - self.boundary_zlo)

                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()
                            
                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Atoms" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if ((len(tokens) <= self.i_atomid) or
                                    (len(tokens) <= self.i_atomtype) or
                                    ((self.i_molid != None) and
                                     (len(tokens) <= self.i_molid))):
                                    raise InputError('Error: The number of columns in the \"Atoms\" section does\n'
                                                     '       not match the atom_style (see column name list above).\n')
                                elif ((len(tokens) != len(self.atom_column_names)) and
                                      (len(tokens) != len(self.atom_column_names) + 3) and
                                      (not self.complained_atom_style_mismatch)):
                                    self.complained_atom_style_mismatch = True
                                    sys.stderr.write('Warning: The number of columns in the \"Atoms\" section does\n'
                                                     '         not match the atom_style (see column name list above).\n')
                                    # this is not a very serious warning.
                                    # self.no_warnings = False <--no need. commenting
                                    # out

                                atomid = Intify(tokens[self.i_atomid])
                                atomtype = Intify(tokens[self.i_atomtype])

                                molid = None
                                if self.i_molid:
                                    molid = Intify(tokens[self.i_molid])

                                self.atomid2type[atomid] = atomtype
                                if self.i_molid:
                                    self.atomid2mol[atomid] = molid

                                if (BelongsToSel(atomid, self.atomid_selection) and
                                    BelongsToSel(atomtype, self.atomtype_selection) and
                                    BelongsToSel(molid, self.molid_selection)):

                                    tokens[self.i_atomid] = '$atom:id' + \
                                        tokens[self.i_atomid]
                                    #tokens[self.i_atomid] = '$atom:'+self.atomids_int2name[atomid]
                                    # fill self.atomtypes_int2str[] with a default name (change later):
                                    #tokens[self.i_atomtype] = '@atom:type'+tokens[self.i_atomtype]
                                    atomtype_name = 'type' + tokens[self.i_atomtype]
                                    tokens[self.i_atomtype] = '@atom:' + atomtype_name

                                    # Interpreting unit-cell counters
                                    # If present, then unit-cell "flags" must be
                                    # added to the x,y,z coordinates.
                                    #
                                    # For more details on unit-cell "flags", see:
                                    # http://lammps.sandia.gov/doc/read_data.html
                                    # "In the data file, atom lines (all lines or
                                    #  none of them) can optionally list 3 trailing
                                    #  integer values (nx,ny,nz), which are used to
                                    #  initialize the atom’s image flags.
                                    #  If nx,ny,nz values are not listed in the
                                    #  data file, LAMMPS initializes them to 0.
                                    #  Note that the image flags are immediately
                                    #  updated if an atom’s coordinates need to
                                    #  wrapped back into the simulation box."

                                    if (len(tokens) == len(self.atom_column_names) + 3):
                                        nx = int(tokens[-3])
                                        ny = int(tokens[-2])
                                        nz = int(tokens[-1])
                                        x = (float(tokens[self.i_x]) +
                                             nx * avec[0] +
                                             ny * bvec[0] +
                                             nz * cvec[0])
                                        y = (float(tokens[self.i_y]) +
                                             nx * avec[1] +
                                             ny * bvec[1] +
                                             nz * cvec[1])
                                        z = (float(tokens[self.i_z]) +
                                             nx * avec[2] +
                                             ny * bvec[2] +
                                             nz * cvec[2])
                                        tokens[self.i_x] = str(x)
                                        tokens[self.i_y] = str(y)
                                        tokens[self.i_z] = str(z)
                                        # Now get rid of them:
                                        del tokens[-3:]

                                    # I can't use self.atomids_int2name or self.atomtypes_int2name yet
                                    # because they probably have not been defined yet.
                                    # (Instead assign these names in a later pass.)

                                    if self.i_molid:
                                        tokens[self.i_molid] = '$mol:m' + \
                                            tokens[self.i_molid]
                                    self.l_data_atoms.append(
                                        (' ' * self.indent) + (' '.join(tokens)))
                                    self.needed_atomids.add(atomid)

                                    self.needed_atomtypes.add(atomtype)
                                    # Not all atom_styles have molids.
                                    # Check for this before adding.
                                    if molid != None:
                                        self.needed_molids.add(molid)

                        for atomtype in self.needed_atomtypes:
                            assert(type(atomtype) is int)
                            if ((self.min_needed_atomtype == None) or
                                (self.min_needed_atomtype > atomtype)):
                                self.min_needed_atomtype = atomtype
                            if ((self.max_needed_atomtype == None) or
                                (self.max_needed_atomtype < atomtype)):
                                self.max_needed_atomtype = atomtype

                        for atomid in self.needed_atomids:
                            assert(type(atomid) is int)
                            if ((self.min_needed_atomid == None) or
                                (self.min_needed_atomid > atomid)):
                                self.min_needed_atomid = atomid
                            if ((self.max_needed_atomid == None) or
                                (self.max_needed_atomid < atomid)):
                                self.max_needed_atomid = atomid
                        for molid in self.needed_molids:
                            assert(type(molid) is int)
                            if ((self.min_needed_molid == None) or
                                (self.min_needed_molid > molid)):
                                self.min_needed_molid = molid
                            if ((self.max_needed_molid == None) or
                                (self.max_needed_molid < molid)):
                                self.max_needed_molid = molid

                        pass

                    elif (line == 'Masses'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            # Read the next line of text but don't skip comments
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  #<- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Masses" section.  Offending line:\n'
                                                     '      "'+line.strip()+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                atomtype = Intify(tokens[0])

                                if BelongsToSel(atomtype, self.atomtype_selection):
                                    tokens[0] = '@atom:type' + tokens[0]
                                    self.l_data_masses.append((' ' * self.indent) +
                                                         (' '.join(tokens)))

                                    # NOTE: We might modify this entry in the
                                    #       self.l_data_masses list later if we
                                    #       learn that the user has chosen a
                                    #       preferred name for this atom type
                                    #       (if self.infer_types_from_comments==True)


                    elif (line == 'Velocities'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Velocities" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                atomid = Intify(tokens[0])
                                atomtype = None
                                if atomid in self.atomid2type:
                                    atomtype = self.atomid2type[atomid]
                                moldid = None
                                if atomid in self.atomid2mol:
                                    molid = self.atomid2mol[atomid]
                                if (BelongsToSel(atomid, self.atomid_selection) and
                                    BelongsToSel(atomtype,
                                                 self.atomtype_selection) and
                                    BelongsToSel(molid,
                                                 self.molid_selection)):

                                    tokens[0] = '$atom:id' + tokens[0]
                                    #tokens[0] = '$atom:'+self.atomids_int2name[atomid]
                                    # NOTE:I can't use "self.atomids_int2name" yet because
                                    #     they probably have not been defined yet.
                                    # (Instead assign these names in a later pass.)
                                    self.l_data_velocities.append(
                                        (' ' * self.indent) + (' '.join(tokens)))

                    # non-point-like-particles:
                    elif (line == 'Ellipsoids'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Ellipsoids" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                atomid = Intify(tokens[0])
                                atomtype = None
                                if atomid in self.atomid2type:
                                    atomtype = self.atomid2type[atomid]
                                moldid = None
                                if atomid in self.atomid2mol:
                                    molid = self.atomid2mol[atomid]
                                if (BelongsToSel(atomid,
                                                 self.atomid_selection) and
                                    BelongsToSel(atomtype,
                                                 self.atomtype_selection) and
                                    BelongsToSel(molid, self.molid_selection)):

                                    tokens[0] = '$atom:id' + tokens[0]
                                    #tokens[0] = '$atom:'+self.atomids_int2name[atomid]
                                    # NOTE:I can't use "self.atomids_int2name" yet because
                                    #     they probably have not been defined yet.
                                    # (Instead assign these names in a later pass.)
                                    self.l_data_ellipsoids.append(
                                        (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Lines'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Lines" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                atomid = Intify(tokens[0])
                                atomtype = None
                                if atomid in self.atomid2type:
                                    atomtype = self.atomid2type[atomid]
                                moldid = None
                                if atomid in self.atomid2mol:
                                    molid = self.atomid2mol[atomid]
                                if (BelongsToSel(atomid, self.atomid_selection) and
                                        BelongsToSel(atomtype, self.atomtype_selection) and
                                        BelongsToSel(molid, self.molid_selection)):
                                    tokens[0] = '$atom:id' + tokens[0]
                                    #tokens[0] = '$atom:'+self.atomids_int2name[atomid]
                                    # NOTE:I can't use "self.atomids_int2name" yet because
                                    #     they probably have not been defined yet.
                                    # (Instead assign these names in a later pass.)
                                    self.l_data_lines.append(
                                        (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Triangles'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Triangles" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                atomid = Intify(tokens[0])
                                atomtype = None
                                if atomid in self.atomid2type:
                                    atomtype = self.atomid2type[atomid]
                                moldid = None
                                if atomid in self.atomid2mol:
                                    molid = self.atomid2mol[atomid]
                                if (BelongsToSel(atomid, self.atomid_selection) and
                                        BelongsToSel(atomtype, self.atomtype_selection) and
                                        BelongsToSel(molid, self.molid_selection)):
                                    tokens[0] = '$atom:id' + tokens[0]
                                    #tokens[0] = '$atom:'+self.atomids_int2name[atomid]
                                    # NOTE:I can't use "self.atomids_int2name" yet because
                                    #     they probably have not been defined yet.
                                    # (Instead assign these names in a later pass.)
                                    self.l_data_triangles.append(
                                        (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Bonds'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Bonds" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 4):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Bonds section:\n'
                                                     '       \"' + line + '\"\n')

                                #tokens[0] = '$bond:id'+tokens[0]
                                #tokens[1] = '@bond:type'+tokens[1]
                                atomids = [None, None]
                                atomtypes = [None, None]
                                molids = [None, None]
                                in_selections = True
                                some_in_selection = False
                                for n in range(0, 2):
                                    atomids[n] = Intify(tokens[2 + n])
                                    if atomids[n] in self.atomid2type:
                                        atomtypes[n] = self.atomid2type[atomids[n]]
                                    if atomids[n] in self.atomid2mol:
                                        molids[n] = self.atomid2mol[atomids[n]]
                                    if (BelongsToSel(atomids[n], self.atomid_selection) and
                                            BelongsToSel(atomtypes[n], self.atomtype_selection) and
                                            BelongsToSel(molids[n], self.molid_selection)):
                                        some_in_selection = True
                                    else:
                                        in_selections = False
                                if in_selections:
                                    self.l_data_bonds.append(
                                        (' ' * self.indent) + (' '.join(tokens)))
                                elif some_in_selection:
                                    sys.stderr.write(
                                        'WARNING: SELECTION BREAKS BONDS\n')
                                    sys.stderr.write(
                                        '         (between atom ids: ')

                                    for n in range(0, 2):
                                        sys.stderr.write(str(atomids[n]) + ' ')
                                    sys.stderr.write(')\n'
                                                     '         The atoms you selected are bonded\n'
                                                     '         to other atoms you didn\'t select.\n'
                                                     '         Are you sure you selected the correct atoms?\n')
                                    self.no_warnings = False

                    elif (line == 'Angles'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Angles" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')


                                if (len(tokens) < 5):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Angles section:\n'
                                                     '       \"' + line + '\"\n')
                                #tokens[0] = '$angle:id'+tokens[0]
                                #tokens[1] = '@angle:type'+tokens[1]
                                atomids = [None, None, None]
                                atomtypes = [None, None, None]
                                molids = [None, None, None]
                                in_selections = True
                                some_in_selection = False
                                for n in range(0, 3):
                                    atomids[n] = Intify(tokens[2 + n])
                                    if atomids[n] in self.atomid2type:
                                        atomtypes[n] = self.atomid2type[atomids[n]]
                                    if atomids[n] in self.atomid2mol:
                                        molids[n] = self.atomid2mol[atomids[n]]
                                    if (BelongsToSel(atomids[n], self.atomid_selection) and
                                            BelongsToSel(atomtypes[n], self.atomtype_selection) and
                                            BelongsToSel(molids[n], self.molid_selection)):
                                        some_in_selection = True
                                    else:
                                        in_selections = False
                                if in_selections:
                                    if not self.ignore_angles_dihedrals_impropers:
                                        self.l_data_angles.append(
                                            (' ' * self.indent) + (' '.join(tokens)))
                                elif some_in_selection:
                                    sys.stderr.write(
                                        'WARNING: SELECTION BREAKS ANGLE INTERACTION\n')
                                    sys.stderr.write(
                                        '         (between atom ids: ')
                                    for n in range(0, 3):
                                        sys.stderr.write(str(atomids[n]) + ' ')
                                    sys.stderr.write(')\n'
                                                     '         The atoms you selected participate in 3-body \"Angle\"\n'
                                                     '         interactions with other atoms you didn\'t select.\n'
                                                     '         (They will be ignored.)\n'
                                                     '         Are you sure you selected the correct atoms?\n')
                                    self.no_warnings = False

                    elif (line == 'Dihedrals'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Dihedrals" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 6):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Dihedrals section:\n'
                                                     '       \"' + line + '\"\n')
                                #tokens[0] = '$dihedral:id'+tokens[0]
                                #tokens[1] = '@dihedral:type'+tokens[1]
                                atomids = [None, None, None, None]
                                atomtypes = [None, None, None, None]
                                molids = [None, None, None, None]
                                in_selections = True
                                some_in_selection = False
                                for n in range(0, 4):
                                    atomids[n] = Intify(tokens[2 + n])
                                    if atomids[n] in self.atomid2type:
                                        atomtypes[n] = self.atomid2type[atomids[n]]
                                    if atomids[n] in self.atomid2mol:
                                        molids[n] = self.atomid2mol[atomids[n]]
                                    if (BelongsToSel(atomids[n], self.atomid_selection) and
                                            BelongsToSel(atomtypes[n], self.atomtype_selection) and
                                            BelongsToSel(molids[n], self.molid_selection)):
                                        some_in_selection = True
                                    else:
                                        in_selections = False
                                if in_selections:
                                    if not self.ignore_angles_dihedrals_impropers:
                                        self.l_data_dihedrals.append(
                                            (' ' * self.indent) + (' '.join(tokens)))
                                elif some_in_selection:
                                    sys.stderr.write(
                                        'WARNING: SELECTION BREAKS DIHEDRAL INTERACTION\n')
                                    sys.stderr.write(
                                        '         (between atom ids: ')
                                    for n in range(0, 4):
                                        sys.stderr.write(str(atomids[n]) + ' ')
                                    sys.stderr.write(')\n'
                                                     '         The atoms you selected participate in 4-body \"Dihedral\"\n'
                                                     '         interactions with other atoms you didn\'t select.\n'
                                                     '         (They will be ignored.)\n'
                                                     '         Are you sure you selected the correct atoms?\n')
                                    self.no_warnings = False


                    elif (line == 'Impropers'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Impropers" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 6):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Impropers section:\n'
                                                     '       \"' + line + '\"\n')
                                #tokens[0] = '$improper:id'+tokens[0]
                                #tokens[1] = '@improper:type'+tokens[1]
                                atomids = [None, None, None, None]
                                atomtypes = [None, None, None, None]
                                molids = [None, None, None, None]
                                in_selections = True
                                some_in_selection = False
                                for n in range(0, 4):
                                    atomids[n] = Intify(tokens[2 + n])
                                    if atomids[n] in self.atomid2type:
                                        atomtypes[n] = self.atomid2type[atomids[n]]
                                    if atomids[n] in self.atomid2mol:
                                        molids[n] = self.atomid2mol[atomids[n]]
                                    if (BelongsToSel(atomids[n], self.atomid_selection) and
                                            BelongsToSel(atomtypes[n], self.atomtype_selection) and
                                            BelongsToSel(molids[n], self.molid_selection)):
                                        some_in_selection = True
                                    else:
                                        in_selections = False
                                if in_selections:
                                    if not self.ignore_angles_dihedrals_impropers:
                                        self.l_data_impropers.append(
                                            (' ' * self.indent) + (' '.join(tokens) + '\n'))
                                elif some_in_selection:
                                    sys.stderr.write(
                                        'WARNING: SELECTION BREAKS IMPROPER INTERACTION\n')
                                    sys.stderr.write(
                                        '         (between atom ids: ')
                                    for n in range(0, 4):
                                        sys.stderr.write(str(atomids[n]) + ' ')
                                    sys.stderr.write(')\n'
                                                     '         The atoms you selected participate in 4-body \"Improper\"\n'
                                                     '         interactions with other atoms you didn\'t select.\n'
                                                     '         (They will be ignored.)\n'
                                                     '         Are you sure you selected the correct atoms?\n')
                                    self.no_warnings = False

                    elif (line == 'CMAP'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "CMAP" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 7):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in CMAP section:\n'
                                                     '       \"' + line + '\"\n')
                                #tokens[0] = '$cmap:id'+tokens[0]
                                #tokens[1] = '@cmap:type'+tokens[1]
                                atomids = [None, None, None, None, None]
                                atomtypes = [None, None, None, None, None]
                                molids = [None, None, None, None, None]
                                in_selections = True
                                some_in_selection = False
                                for n in range(0, 5):
                                    atomids[n] = Intify(tokens[2 + n])
                                    if atomids[n] in self.atomid2type:
                                        atomtypes[n] = self.atomid2type[atomids[n]]
                                    if atomids[n] in self.atomid2mol:
                                        molids[n] = self.atomid2mol[atomids[n]]
                                    if (BelongsToSel(atomids[n], self.atomid_selection) and
                                            BelongsToSel(atomtypes[n], self.atomtype_selection) and
                                            BelongsToSel(molids[n], self.molid_selection)):
                                        some_in_selection = True
                                    else:
                                        in_selections = False
                                if in_selections:
                                    if not self.ignore_angles_dihedrals_impropers:
                                        self.l_data_cmap.append(
                                            (' ' * self.indent) + (' '.join(tokens)))
                                elif some_in_selection:
                                    sys.stderr.write(
                                        'WARNING: SELECTION BREAKS CMAP INTERACTION\n')
                                    sys.stderr.write(
                                        '         (between atom ids: ')
                                    for n in range(0, 4):
                                        sys.stderr.write(str(atomids[n]) + ' ')
                                    sys.stderr.write(')\n'
                                                     '         The atoms you selected participate in 5-body \"CMAP\"\n'
                                                     '         interactions with other atoms you didn\'t select.\n'
                                                     '         (They will be ignored.)\n'
                                                     '         Are you sure you selected the correct atoms?\n')
                                    self.no_warnings = False

                    elif (line == 'Bond Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Bond Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@bond:type'+tokens[0]
                                self.l_data_bond_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Angle Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:
                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Angle Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@angle:type'+tokens[0]
                                self.l_data_angle_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Dihedral Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Dihedral Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_dihedral_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Improper Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Improper Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@improper:type'+tokens[0]
                                self.l_data_improper_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Pair Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        self.some_pair_coeffs_read = True
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Pair Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 2):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Pair Coeffs section:\n'
                                                     '       \"' + line + '\"\n')
                                atomtype_i_str = tokens[0]
                                if '*' in atomtype_i_str:
                                    raise InputError('PROBLEM near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '         As of 2017-10, moltemplate forbids use of the "\*\" wildcard\n'
                                                     '         character in the \"Pair Coeffs\" section.\n')
                                else:
                                    i = int(atomtype_i_str)
                                    if ((not i) or
                                            BelongsToSel(i, self.atomtype_selection)):
                                        ##i_str = '@atom:type' + str(i)
                                        #i_str = '@atom:' + self.atomtypes_int2name[i]
                                        i_str = str(i)
                                        tokens[0] = i_str
                                        self.l_data_pair_coeffs.append(
                                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'PairIJ Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        self.some_pair_coeffs_read = True
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "PairIJCoeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                if (len(tokens) < 2):
                                    raise InputError('Error: near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '       Nonsensical line in Pair Coeffs section:\n'
                                                     '       \"' + line + '\"\n')
                                atomtype_i_str = tokens[0]
                                atomtype_j_str = tokens[1]
                                if (('*' in atomtype_i_str) or ('*' in atomtype_j_str)):
                                    raise InputError('PROBLEM near or before ' +
                                                     ErrorLeader(lex.infile,
                                                                 lex.lineno) + '\n'
                                                     '         As of 2017-10, moltemplate forbids use of the "\*\" wildcard\n'
                                                     '         character in the \"PairIJ Coeffs\" section.\n')
                                else:
                                    i = int(atomtype_i_str)
                                    j = int(atomtype_j_str)
                                    if (((not i) or BelongsToSel(i, self.atomtype_selection)) and
                                            ((not j) or BelongsToSel(j, self.atomtype_selection))):
                                        ##i_str = '@atom:type' + str(i)
                                        ##j_str = '@atom:type' + str(j)
                                        #i_str = '@atom:' + self.atomtypes_int2name[i]
                                        #j_str = '@atom:' + self.atomtypes_int2name[j]
                                        i_str = str(i)
                                        j_str = str(j)
                                        tokens[0] = i_str
                                        tokens[1] = j_str
                                        self.l_data_pairij_coeffs.append(
                                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'pair_coeff'):
                        self.some_pair_coeffs_read = True
                        if (len(tokens) < 3):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical pair_coeff command:\n'
                                             '       \"' + line + '\"\n')
                        self.l_in_pair_coeffs.append(' ' * self.indent + line)

                    elif (tokens[0] == 'mass'):
                        self.some_pair_coeffs_read = True
                        if (len(tokens) < 3):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical \"mass\" command:\n'
                                             '       \"' + line + '\"\n')
                        self.l_in_masses.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'bond_coeff'):
                        if (len(tokens) < 2):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical bond_coeff command:\n'
                                             '       \"' + line + '\"\n')
                        #tokens[1] = '@bond:type'+tokens[1]
                        self.l_in_bond_coeffs.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'angle_coeff'):
                        if (len(tokens) < 2):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical angle_coeff command:\n'
                                             '       \"' + line + '\"\n')
                        #tokens[1] = '@angle:type'+tokens[1]
                        self.l_in_angle_coeffs.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'dihedral_coeff'):
                        if (len(tokens) < 2):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical dihedral_coeff command:\n'
                                             '       \"' + line + '\"\n')
                        #tokens[1] = '@dihedral:type'+tokens[1]
                        self.l_in_dihedral_coeffs.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'improper_coeff'):
                        if (len(tokens) < 2):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical improper_coeff command:\n'
                                             '       \"' + line + '\"\n')
                        #tokens[1] = '@improper:type'+tokens[1]
                        self.l_in_improper_coeffs.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    # -- class2 force fields --
                    elif (line == 'BondBond Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "BondBond Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@angle:type'+tokens[0]
                                self.l_data_bondbond_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'BondAngle Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "BondAngle Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@angle:type'+tokens[0]
                                self.l_data_bondangle_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'MiddleBondTorsion Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "MiddleBondTorsion Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_middlebondtorsion_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'EndBondTorsion Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "EndBondTorsion Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_endbondtorsion_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'AngleTorsion Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "AngleTorsion Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_angletorsion_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'AngleAngleTorsion Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "AngleAngleTorsion Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_angleangletorsion_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'BondBond13 Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "BondBond13 Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@dihedral:type'+tokens[0]
                                self.l_data_bondbond13_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'AngleAngle Coeffs'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "AngleAngle Coeffs" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                #tokens[0] = '@improper:type'+tokens[0]
                                self.l_data_angleangle_coeffs.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Angles By Type'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Angles By Type" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                tokens[0] = '@angle:type' + tokens[0]
                                self.l_data_angles_by_type.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Dihedrals By Type'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Dihedrals By Type" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                tokens[0] = '@dihedral:type' + tokens[0]
                                self.l_data_dihedrals_by_type.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    elif (line == 'Impropers By Type'):
                        sys.stderr.write('  reading \"' + line + '\"\n')
                        while lex:
                            line_orig = lex.ReadLine()
                            line = line_orig
                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.strip()

                            if line in data_file_header_names:
                                lex.push_raw_text(line_orig+'\n')  # <- Save line for later
                                break

                            tokens = line.split()

                            if len(tokens) > 0:

                                if not tokens[0].strip().isdigit():
                                    raise InputError('Error near the "Impropers By Type" section.  Offending line:\n'
                                                     '      "'+line+'"\n'
                                                     '      This line does not begin with an integer, and it does not appear to be\n'
                                                     '      a LAMMPS data file section name. If this is a valid data section name,\n'
                                                     '      please let the author know.  (jewett.aij -at- gmail.com)\n')

                                tokens[0] = '@improper:type' + tokens[0]
                                self.l_data_impropers_by_type.append(
                                    (' ' * self.indent) + (' '.join(tokens)))

                    # Figure out the size of the simulation box boundary:
                    elif ((len(tokens) == 4) and
                          (tokens[2] == 'xlo') and
                          (tokens[3] == 'xhi') and
                          IsNumber(tokens[0]) and
                          IsNumber(tokens[1])):
                        self.boundary_xlo = float(tokens[0])
                        self.boundary_xhi = float(tokens[1])

                    elif ((len(tokens) == 4) and
                          (tokens[2] == 'ylo') and
                          (tokens[3] == 'yhi') and
                          IsNumber(tokens[0]) and
                          IsNumber(tokens[1])):
                        self.boundary_ylo = float(tokens[0])
                        self.boundary_yhi = float(tokens[1])

                    elif ((len(tokens) == 4) and
                          (tokens[2] == 'zlo') and
                          (tokens[3] == 'zhi') and
                          IsNumber(tokens[0]) and
                          IsNumber(tokens[1])):
                        self.boundary_zlo = float(tokens[0])
                        self.boundary_zhi = float(tokens[1])

                    elif ((len(tokens) == 6) and
                          (tokens[3] == 'xy') and
                          (tokens[4] == 'xz') and
                          (tokens[5] == 'yz') and
                          IsNumber(tokens[0]) and
                          IsNumber(tokens[1]) and
                          IsNumber(tokens[2])):
                        self.boundary_xy = float(tokens[0])
                        self.boundary_xz = float(tokens[1])
                        self.boundary_yz = float(tokens[2])

                    elif (tokens[0] == 'group'):
                        if (len(tokens) < 3):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical group command:\n'
                                             '       \"' + line.strip() + '\"\n')
                        self.l_in_group.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif (tokens[0] == 'set'):
                        if (len(tokens) < 3):
                            raise InputError('Error: near or before ' +
                                             ErrorLeader(lex.infile,
                                                         lex.lineno) + '\n'
                                             '       Nonsensical set command:\n'
                                             '       \"' + line.strip() + '\"\n')
                        self.l_in_set.append(
                            (' ' * self.indent) + (' '.join(tokens)))

                    elif ((tokens[0] == 'fix') and (len(tokens) >= 4)):
                        if (tokens[3].find('rigid') == 0):
                            if (len(tokens) < 6):
                                raise InputError('Error: near or before ' +
                                                 ErrorLeader(lex.infile,
                                                             lex.lineno) + '\n'
                                                 '       Nonsensical ' +
                                                 tokens[0] + ' ' +
                                                 tokens[3] + ' command:\n'
                                                 '       \"' + line.strip() + '\"\n')
                            self.l_in_fix_rigid.append(
                                (' ' * self.indent) + (' '.join(tokens)))
                        elif (tokens[3].find('shake') == 0):
                            if (len(tokens) < 7):
                                raise InputError('Error: near or before ' +
                                                 ErrorLeader(lex.infile,
                                                             lex.lineno) + '\n'
                                                 '       Nonsensical ' +
                                                 tokens[0] + ' ' +
                                                 tokens[3] + ' command:\n'
                                                 '       \"' + line.strip() + '\"\n')
                            self.l_in_fix_shake.append(
                                (' ' * self.indent) + (' '.join(tokens)))
                        elif (tokens[3].find('poems') == 0):
                            if (len(tokens) < 4):
                                raise InputError('Error: near or before ' +
                                                 ErrorLeader(lex.infile,
                                                             lex.lineno) + '\n'
                                                 '       Nonsensical ' +
                                                 tokens[0] + ' ' +
                                                 tokens[3] + ' command:\n'
                                                 '       \"' + line.strip() + '\"\n')
                            self.l_in_fix_poems.append(
                                (' ' * self.indent) + (' '.join(tokens)))
                        elif (tokens[3].find('qeq') == 0):
                            if (len(tokens) < 8):
                                raise InputError('Error: near or before ' +
                                                 ErrorLeader(lex.infile,
                                                             lex.lineno) + '\n'
                                                 '       Nonsensical ' +
                                                 tokens[0] + ' ' +
                                                 tokens[3] + ' command:\n'
                                                 '       \"' + line.strip() + '\"\n')
                            self.l_in_fix_qeq.append(
                                (' ' * self.indent) + (' '.join(tokens)))
                        elif (tokens[3].find('qmmm') == 0):
                            if (len(tokens) < 8):
                                raise InputError('Error: near or before ' +
                                                 ErrorLeader(lex.infile,
                                                             lex.lineno) + '\n'
                                                 '       Nonsensical ' +
                                                 tokens[0] + ' ' +
                                                 tokens[3] + ' command:\n'
                                                 '       \"' + line.strip() + '\"\n')
                            self.l_in_fix_qmmm.append(
                                (' ' * self.indent) + (' '.join(tokens)))
                        elif (tokens[3].find('restrain') == 0):
                            sys.stderr.write('WARNING: fix \"' + tokens[3] + '\" commands are NOT understood by ' + g_program_name + '.\n'
                                       '  If you need restraints, add them to your final .LT file (eg. \"system.lt\"),\n'
                                       '  (And be sure to use unique (full, long) moltemplate names for each $atom:.)\n'
                                       '  Ignoring line \"' + line.strip() + '\"\n')


                    elif (self.infer_types_from_comments and
                          (len(tokens) == 3) and (tokens[2] == 'types') and
                          ((tokens[1] == 'atom') or
                           (tokens[1] == 'bond') or
                           (tokens[1] == 'angle') or
                           (tokens[1] == 'dihedral') or
                           (tokens[1] == 'improper'))):


                        #####################################################
                        # Consider comments following these commands:
                        #"N atom types"
                        #"N bond types"
                        #"N angle types"
                        #"N dihedral types"
                        #"N improper types"
                        # Lines containing nothing but comments after
                        # these lines are ASSUMED to contain names
                        # for these bonded types.
                        # So instead of generating bond types with names like:
                        # @bond:type1, @bond:type2, you will get names like:
                        # @bond:CA_HC1, @{bond:N CA}, ...
                        #####################################################

                        interaction_type = tokens[1]
                        #interaction_type='bond','angle','dihedral','improper', (or 'atom')

                        ntypes = 0
                        encountered_types = set([])

                        sys.stderr.write('  reading \"' + line.strip() + '\"\n')
                        while lex:
                            # Read the next line of text but don't skip comments
                            #comment_char_backup = lex.commenters
                            lex.commenters = ''
                            line = lex.ReadLine()
                            #lex.commenters = comment_char_backup

                            comment_text = ''
                            ic = line.find('#')
                            if ic != -1:
                                comment_text = line[ic + 1:].strip()
                                line = line[:ic]
                            line = line.rstrip()

                            #ltokens = SplitQuotedString(line.strip())
                            ltokens = line.strip().split()

                            if ((line.strip() in data_file_header_names) or
                                (len(ltokens) >= 2) and
                                ltokens[0].isdigit() and
                                ((len(ltokens) == 3) and
                                 (ltokens[1] == 'atom') and
                                 (ltokens[2] == 'types')) or
                                ((len(ltokens) == 3) and
                                 (ltokens[1] == 'bond') and
                                 (ltokens[2] == 'types')) or
                                ((len(ltokens) == 3) and
                                 (ltokens[1] == 'angle') and
                                 (ltokens[2] == 'types')) or
                                ((len(ltokens) == 3) and
                                 (ltokens[1] == 'dihedral') and
                                 (ltokens[2] == 'types')) or
                                ((len(ltokens) == 3) and
                                 (ltokens[1] == 'improper') and
                                 (ltokens[2] == 'types')) or
                                ((len(ltokens) == 2) and
                                 (ltokens[1] == 'crossterms')) or
                                ((len(ltokens) == 4) and
                                 (ltokens[2] == 'xlo') and
                                 (ltokens[3] == 'xhi')) or
                                ((len(ltokens) == 4) and
                                 (ltokens[2] == 'ylo') and
                                 (ltokens[3] == 'yhi')) or
                                ((len(ltokens) == 4) and
                                 (ltokens[2] == 'zlo') and
                                 (ltokens[3] == 'zhi')) or
                                ((len(ltokens) == 6) and
                                 (ltokens[3] == 'xy') and
                                 (ltokens[4] == 'xz') and
                                 (ltokens[5] == 'yz'))):
                                lex.push_raw_text(line+'\n')  #<- Save line for later
                                break

                            #ctokens = SplitQuotedString(comment_text.strip())
                            #if len(ctokens) == 1:
                            #    type_int = ntypes+1
                            #    type_str = RemoveOuterQuotes(ctokens[0],'\'\"')
                            #elif len(ctokens) >= 2:
                            #    type_int = Intify(ctokens[0])
                            #    type_str = RemoveOuterQuotes(ctokens[1],'\'\"')
                            #else:
                            #    continue

                            if len(ltokens) > 0:
                                continue

                            comment_text = comment_text.strip()
                            ctokens = comment_text.split()
                            if len(ctokens) == 1:
                                type_int = ntypes+1
                                type_str = ctokens[0]
                            elif len(ctokens) >= 2:
                                type_int = Intify(ctokens[0])
                                #type_str = ctokens[1]

                                #iws=comment_text.find(' ') <--only works on ' '
                                iws = 0 #<--detect first character of whitespace
                                while iws < len(comment_text):
                                    if comment_text[iws].isspace():
                                        break
                                    iws += 1
                                type_str = comment_text[iws+1:].strip()
                            else:
                                continue

                            if (interaction_type == 'atom'):
                                self.atomtypes_int2name[type_int] = type_str
                                self.atomtypes_name2int[type_str] = type_int
                                ntypes += 1
                            if (interaction_type == 'bond'):
                                self.bondtypes_int2name[type_int]=type_str
                                ntypes += 1
                            elif (interaction_type == 'angle'):
                                self.angletypes_int2name[type_int]=type_str
                                ntypes += 1
                            elif (interaction_type == 'dihedral'):
                                self.dihtypes_int2name[type_int]=type_str
                                ntypes += 1
                            elif (interaction_type == 'improper'):
                                self.imptypes_int2name[type_int]=type_str
                                ntypes += 1
                            else:
                                continue
                            if (self.forbid_type_name_duplicates and
                                (type_str in encountered_types)):
                                raise InputError('Error: Duplicate name in comments following: \"' + interaction_type + ' types\" section\n'
                                                 '       (By default ' + g_program_name +
                                                 ' attempts to infer '+interaction_type+' type names from\n'
                                                 '       comments following the \"'+interaction_type+'\" section of your data file.)\n'
                                                 '       This error message occurs if two different '+interaction_type+' types\n'
                                                 '       have the same name (in the comments section)\n'
                                                 '       You can avoid this error by adding the \"-ignore-duplicates\" argument.\n')


                            # Print out an optional debug string:
                            #sys.stderr.write(interaction_type+'_type_name['+
                            #                 str(type_int)+'] = \"'+
                            #                 type_str+'\"\n');
                                            

                    else:
                        sys.stderr.write('  Ignoring line \"' +
                                         line.strip() + '\"\n')


    def PostProcess1(self):
        """
        PostProcess1:
        Now do a second-pass throught the "self.l_data_atoms" section, and
        finish dealing with "self.infer_types_from_comments".

        Determine atomid names from atomtype names.

        If an atom type appears multiple times (as is usually the case),
        then set the atomid names to the atomtype name with an integer suffix.
        (For example, the atoms of type "CA" will be named
         "CA_1", "CA_2", "CA_3", ...)

        If any atom types only appear once,
        set the atomid name = atomtype name.
        """
        
        sys.stderr.write('\n\n')

        sys.stderr.write('  processing \"Atoms\" section (\n')
        sys.stderr.write('    postprocess1,\n')


        # set the atomid name to the atomtype name with an integer suffix.
        # (For example, the atoms of type "CA" will be named
        #  "CA_1", "CA_2", "CA_3", ...)

        for i in range(0, len(self.l_data_atoms)):
            tokens = self.l_data_atoms[i].split()
            atomid = tokens[self.i_atomid]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)

            simple_atom_type_names = True
            atomtype = tokens[self.i_atomtype]
            # remove the "@atom:" prefix (we will put it back later)
            if atomtype.find('@atom:') == 0:
                atomtype = Intify(atomtype[6:])
                # convert to an integer
            if not (atomtype in self.atomtypes_int2name):
                self.atomtypes_int2name[atomtype] = 'type' + str(atomtype)
            if self.infer_types_from_comments:
                atomtype_name = self.atomtypes_int2name[atomtype]
                if atomtype_name.find('/') != -1:
                    simple_atom_type_names = False
            if (self.infer_types_from_comments and simple_atom_type_names):
                if atomtype in self.atomids_by_type:
                    self.l_atomids = self.atomids_by_type[atomtype]
                    prev_count = len(self.l_atomids)
                    # lookup the most recently added atom of this type:
                    #prev_atomid_name = self.l_atomids[-1]
                    #ic = prev_atomid_name.rfind('_')
                    #prev_count = int(prev_atomid_name[ic+1:])
                    atomid_name = atomtype_name + '_' + str(prev_count + 1)
                else:
                    atomid_name = atomtype_name  + '_1'
                self.atomids_int2name[atomid] = atomid_name
                #self.atomids_name2str[atomid_name] = atomid
            else:
                self.atomids_int2name[atomid] = 'id' + str(atomid)
            self.atomids_by_type[atomtype].append(atomid)

        # If any atom types only appear once, simplify atomid names
        # by setting the atomid name = atomtype name.

        for i in range(0, len(self.l_data_atoms)):
            tokens = self.l_data_atoms[i].split()

            # remove the "@atom:" prefix (we will put it back later)
            atomtype = tokens[self.i_atomtype]
            if atomtype.find('@atom:') == 0:
                atomtype = atomtype[6:]
            atomtype = Intify(atomtype)
            if self.infer_types_from_comments:
                if len(self.atomids_by_type[atomtype]) == 1:
                    atomid = tokens[self.i_atomid]
                    if atomid.find('$atom:') == 0:
                        atomid = atomid[6:]
                    atomid = Intify(atomid)
                    atomtype_name = self.atomtypes_int2name[atomtype]
                    self.atomids_int2name[atomid] = atomtype_name


    def PostProcess2(self):
        """
        PostProcess2:
             Substitute the updated atomid names and atom type names into:
        self.l_data_atoms
        self.l_data_velocities
        self.l_data_ellipsoids
        self.l_data_lines
        self.l_data_triangles
        """

        sys.stderr.write('    postprocess2,\n')
        
        for i in range(0, len(self.l_data_atoms)):
            tokens = self.l_data_atoms[i].split()
            atomid = tokens[self.i_atomid]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)
            atomtype = tokens[self.i_atomtype]
            if atomtype.find('@atom:') == 0:
                atomtype = atomtype[6:]
            atomtype = Intify(atomtype)
            tokens = self.l_data_atoms[i].split()
            atomid_name = self.atomids_int2name[atomid]
            atomtype_name = self.atomtypes_int2name[atomtype]
            if atomid_name.find(' ') != -1:
                tokens[self.i_atomid] = '${atom:' + atomid_name + '}'
            else:
                tokens[self.i_atomid] = '$atom:' + atomid_name
            if atomtype_name.find(' ') != -1:
                tokens[self.i_atomtype] = '@{atom:' +atomtype_name+'}'
            else:
                tokens[self.i_atomtype] = '@atom:' + atomtype_name
            self.l_data_atoms[i] = (' ' * self.indent) + (' '.join(tokens))

        # self.l_data_velocities
        for i in range(0, len(self.l_data_velocities)):
            tokens = self.l_data_velocities[i].split()
            atomid = tokens[0]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)
            atomid_name = self.atomids_int2name[atomid]
            tokens = self.l_data_velocities[i].split()
            if atomid_name.find(' ') != -1:
                tokens[0] = '${atom:' + atomid_name + '}'
            else:
                tokens[0] = '$atom:' + atomid_name
            self.l_data_velocities[i] = (' ' * self.indent) + (' '.join(tokens))

        # self.l_data_ellipsoids
        for i in range(0, len(self.l_data_ellipsoids)):
            tokens = self.l_data_ellipsoids[i].split()
            atomid = tokens[0]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)
            atomid_name = self.atomids_int2name[atomid]
            tokens = self.l_data_ellipsoids[i].split()
            if atomid_name.find(' ') != -1:
                tokens[0] = '${atom:' + atomid_name + '}'
            else:
                tokens[0] = '$atom:' + atomid_name
            self.l_data_ellipsoids[i] = (' ' * self.indent) + (' '.join(tokens))

        # self.l_data_lines
        for i in range(0, len(self.l_data_lines)):
            tokens = self.l_data_lines[i].split()
            atomid = tokens[0]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)
            atomid_name = self.atomids_int2name[atomid]
            tokens = self.l_data_lines[i].split()
            if atomid_name.find(' ') != -1:
                tokens[0] = '${atom:' + atomid_name + '}'
            else:
                tokens[0] = '$atom:' + atomid_name
            self.l_data_lines[i] = (' ' * self.indent) + (' '.join(tokens))

        # self.l_data_triangles
        for i in range(0, len(self.l_data_triangles)):
            tokens = self.l_data_triangles[i].split()
            atomid = tokens[0]
            if atomid.find('$atom:') == 0:
                atomid = atomid[6:]
                # convert to an integer
                atomid = Intify(atomid)
            atomid_name = self.atomids_int2name[atomid]
            tokens = self.l_data_triangles[i].split()
            if atomid_name.find(' ') != -1:
                tokens[0] = '${atom:' + atomid_name + '}'
            else:
                tokens[0] = '$atom:' + atomid_name
            self.l_data_triangles[i] = (' ' * self.indent) + (' '.join(tokens))

        if len(self.l_data_atoms) == 0:
            raise InputError('Error(' + g_program_name + '): You have no atoms in you selection!\n'
                             '\n'
                             '       Either you have chosen a set of atoms, molecules, or atom types which\n'
                             '       does not exist, or there is a problem with (the format of) your\n'
                             '       arguments. Check the documentation and examples.\n')






    def PostProcess3(self):
        """
        PostProcess3
        1) Delete lines from l_data_bonds, l_data_angles, l_data_dihedrals,
           l_data_impropers, and cmaps describing interactions between atoms we
           don't care about (i.e. atoms which were not selected by the user).
        2) For the atoms we do care about, replace their atomid names
           (which in this section of the file are still currently numeric)
           with the updated atomids that we determined earlier.
        3) Now that we know which bonds, angles, dihedrals, and impropers
           and cmaps we want to keep, we can figure out which bond types,
           angle types, dihedral types, and improper types are no longer
           needed, and delete their corresponding entry from the lists:
           l_in_bond_coeffs, l_in_angle_coeffs, l_in_dihedral_coeffs, and
           l_in_improper_coeffs
        4) Do the same thing for l_in_pair_coeffs and l_data_masses.  (Discard
           pair_coeffs and masses for atom types we don't care about.)
        """
           
        sys.stderr.write('    postprocess3\n')

        # --- Now delete items that were not selected from the other lists ---

        # --- MASSES ---

        # delete masses for atom types we no longer care about:
        # also substitute the correct atom type name
        i_line = 0
        while i_line < len(self.l_data_masses):
            line = self.l_data_masses[i_line]
            tokens = line.strip().split()
            atomtypestr = tokens[0]
            if atomtypestr.find('@atom:') == 0:
                atomtypestr = atomtypestr[6:]
            atomtype = Intify(atomtypestr)
            if ((not (atomtype in self.needed_atomtypes)) and
                (not ((len(self.atomtype_selection) > 0) and
                      BelongsToSel(atomtype, self.atomtype_selection)))):
                del self.l_data_masses[i_line]
            else:
                tokens[0]=Stringify(atomtype,
                                    self.atomtypes_int2name,
                                    '@','atom','type')
                self.l_data_masses[i_line] = ((' ' * self.indent) +
                                              (' '.join(tokens)))
                i_line += 1

        # --- PAIR COEFFS ---

        # delete data_pair_coeffs for atom types we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_pair_coeffs):
            line = self.l_data_pair_coeffs[i_line]
            tokens = line.strip().split()
            assert(len(tokens) > 0)
            atomtypestr = tokens[0]
            if atomtypestr.find('@atom:') == 0:
                atomtypestr = atomtypestr[6:]
            atomtype = Intify(atomtypestr)
            if ((not (atomtype in self.needed_atomtypes)) and
                (not ((len(self.atomtype_selection) > 0) and
                      BelongsToSel(atomtype, self.atomtype_selection)))):
                del self.l_data_pair_coeffs[i_line]
            else:
                tokens[0] = Stringify(atomtype,
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                self.l_data_pair_coeffs[i_line] = ((' ' * self.indent) +
                                                   (' '.join(tokens)))
                i_line += 1

        # delete data_pairij_coeffs for atom types we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_pairij_coeffs):
            line = self.l_data_pairij_coeffs[i_line]
            tokens = line.strip().split()
            assert(len(tokens) > 0)

            atomtypestr_I = tokens[0]
            if atomtypestr_I.find('@atom:') == 0:
                atomtypestr_I = atomtypestr[6:]
            atomtype_I = Intify(atomtypestr_I)

            atomtypestr_J = tokens[1]
            if atomtypestr_J.find('@atom:') == 0:
                atomtypestr_J = atomtypestr[6:]
            atomtype_J = Intify(atomtypestr_J)

            if (((not (atomtype_I in self.needed_atomtypes)) and
                 (not ((len(self.atomtype_selection) > 0) and
                       BelongsToSel(atomtype_I, self.atomtype_selection))))
                or
                ((not (atomtype_J in self.needed_atomtypes)) and
                 (not ((len(self.atomtype_selection) > 0) and
                       BelongsToSel(atomtype_J, self.atomtype_selection))))):
                del self.l_data_pairij_coeffs[i_line]
            else:
                tokens[0] = Stringify(atomtype_I,
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                tokens[1] = Stringify(atomtype_J,
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                self.l_data_pairij_coeffs[i_line] = ((' ' * self.indent) +
                                                     (' '.join(tokens)))
                i_line += 1

        # delete in_pair_coeffs for atom we no longer care about:
        i_line = 0
        while i_line < len(self.l_in_pair_coeffs):
            line = self.l_in_pair_coeffs[i_line]
            tokens = line.strip().split()
            atomtype_i_str = tokens[1]
            atomtype_j_str = tokens[2]
            #if (('*' in atomtype_i_str) or
            #    ('*' in atomtype_j_str)):
            #    sys.stderr.write('WARNING: near or before '+ErrorLeader(lex.infile, lex.lineno)+'\n'
            #                     '         pair_coeff command contains a \"*\" character.\n'
            #                     '         Keep in mind that using moltemplate.sh you can manually change the\n'
            #                     '         numbers assigned to each atom type (when using -a or -b).  Make sure\n'
            #                     '         nor to accidentally change the order of atom types in one of these\n'
            #                     '         pair_coeff commands.  For example, commands like\n'
            #                     '            pair_coeff 10*4 20*10 0.15 3.6\n'
            #                     '         can be generated by moltemplate.sh, however\n'
            #                     '         they may be rejected by LAMMPS (because LAMMPS prefers this\n'
            #                     '            pair_coeff 4*10 10*20 0.15 3.6)\n'
            #                     '         Later on, you may want to check to make sure moltemplate.sh\n'
            #                     '         is not doing this.  (Fortunately you never have to worry unless\n'
            #                     '         you are using the -a or -b arguments with moltemplate.sh)\n')

            if ('*' in atomtype_i_str):
                atomtype_i_tokens = atomtype_i_str.split('*')

                if atomtype_i_tokens[0] == '':
                    if (self.min_sel_atomtype and
                        (self.min_sel_atomtype < self.min_needed_atomtype)):
                        i_a = self.min_sel_atomtype
                    else:
                        i_a = self.min_needed_atomtype
                else:
                    i_a = Intify(atomtype_i_tokens[0])

                if atomtype_i_tokens[1] == '':
                    if (self.max_sel_atomtype and
                        (self.max_sel_atomtype > self.max_needed_atomtype)):
                        i_b = self.max_sel_atomtype
                    else:
                        i_b = self.max_needed_atomtype
                else:
                    i_b = Intify(atomtype_i_tokens[1])

            else:
                i_a = i_b = Intify(atomtype_i_str)

            assert((type(i_a) is int) and (type(i_b) is int))

            i_a_final = None
            i_b_final = None

            for i in range(i_a, i_b + 1):
                if ((i in self.needed_atomtypes) or
                    ((self.min_sel_atomtype != None) and (self.min_sel_atomtype <= i))):
                    i_a_final = i
                    break

            for i in reversed(range(i_a, i_b + 1)):
                if ((i in self.needed_atomtypes) or
                    ((self.max_sel_atomtype != None) and (self.max_sel_atomtype >= i))):
                    i_b_final = i
                    break

            # if i_a_final and i_b_final:
            #    if i_a_final == i_b_final:
            #        i_str = '@atom:type'+str(i_a_final)
            #        tokens[1] = i_str
            #    else:
            #        i_str = '@{atom:type'+str(i_a_final)+'}*@{atom:type'+str(i_b_final)+'}'

            if ('*' in atomtype_j_str):
                atomtype_j_tokens = atomtype_j_str.split('*')

                if atomtype_j_tokens[0] == '':
                    if (self.min_sel_atomtype and
                        (self.min_sel_atomtype < self.min_needed_atomtype)):
                        j_a = self.min_sel_atomtype
                    else:
                        j_a = self.min_needed_atomtype
                else:
                    j_a = Intify(atomtype_j_tokens[0])

                if atomtype_j_tokens[1] == '':
                    if (self.max_sel_atomtype and
                        (self.max_sel_atomtype > self.max_needed_atomtype)):
                        j_b = self.max_sel_atomtype
                    else:
                        j_b = self.max_needed_atomtype
                else:
                    j_b = Intify(atomtype_j_tokens[1])

            else:
                j_a = j_b = Intify(atomtype_j_str)

            j_a_final = None
            j_b_final = None
            for j in range(j_a, j_b + 1):
                if ((j in self.needed_atomtypes) or
                    ((self.min_sel_atomtype != None) and (self.min_sel_atomtype <= j))):
                    j_a_final = j
                    break
            for j in reversed(range(j_a, j_b + 1)):
                if ((j in self.needed_atomtypes) or
                    ((self.max_sel_atomtype != None) and (self.max_sel_atomtype >= j))):
                    j_b_final = j
                    break

            # if j_a_final and j_b_final:
            #    if j_a_final == j_b_final:
            #        j_str = '@atom:type'+str(j_a_final)
            #        tokens[1] = j_str
            #    else:
            #        j_str = '@{atom:type'+str(j_a_final)+'}*@{atom:type'+str(j_b_final)+'}'

            if not (i_a_final and i_b_final and j_a_final and j_b_final):
                del self.l_in_pair_coeffs[i_line]
            elif (('*' in atomtype_i_str) or ('*' in atomtype_j_str)):
                del self.l_in_pair_coeffs[i_line]
                for i in range(i_a_final, i_b_final + 1):
                    for j in range(j_a_final, j_b_final + 1):
                        if j >= i:
                            tokens[1] = Stringify(i,
                                                  self.atomtypes_int2name,
                                                  '@','atom','type')
                            tokens[2] = Stringify(j,
                                                  self.atomtypes_int2name,
                                                  '@','atom','type')
                            self.l_in_pair_coeffs.insert(i_line,
                                                         (' ' * self.indent) + (' '.join(tokens)))
                            i_line += 1
            else:
                tokens[1] = Stringify(int(tokens[1]),
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                tokens[2] = Stringify(int(tokens[2]),
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                self.l_in_pair_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # delete mass commands for atom types we no longer care about:
        i_line = 0
        while i_line < len(self.l_in_masses):
            line = self.l_in_masses[i_line]
            tokens = line.strip().split()
            atomtype_i_str = tokens[1]
            # if (('*' in atomtype_i_str) or
            #    ('*' in atomtype_j_str)):
            #    sys.stderr.write('WARNING: near or before '+ErrorLeader(lex.infile, lex.lineno)+'\n'
            #                     '         pair_coeff command contains a \"*\" character.\n'
            #                     '         Keep in mind that using moltemplate.sh you can manually change the\n'
            #                     '         numbers assigned to each atom type (when using -a or -b).  Make sure\n'
            #                     '         nor to accidentally change the order of atom types in one of these\n'
            #                     '         pair_coeff commands.  For example, commands like\n'
            #                     '            pair_coeff 10*4 20*10 0.15 3.6\n'
            #                     '         can be generated by moltemplate.sh, however\n'
            #                     '         they may be rejected by LAMMPS (because LAMMPS prefers this\n'
            #                     '            pair_coeff 4*10 10*20 0.15 3.6)\n'
            #                     '         Later on, you may want to check to make sure moltemplate.sh\n'
            #                     '         is not doing this.  (Fortunately you never have to worry unless\n'
            #                     '         you are using the -a or -b arguments with moltemplate.sh)\n')

            if ('*' in atomtype_i_str):
                atomtype_i_tokens = atomtype_i_str.split('*')

                if atomtype_i_tokens[0] == '':
                    if (self.min_sel_atomtype and
                        (self.min_sel_atomtype < self.min_needed_atomtype)):
                        i_a = self.min_sel_atomtype
                    else:
                        i_a = self.min_needed_atomtype
                else:
                    i_a = Intify(atomtype_i_tokens[0])

                if atomtype_i_tokens[1] == '':
                    if (self.max_sel_atomtype and
                        (self.max_sel_atomtype > self.max_needed_atomtype)):
                        i_b = self.max_sel_atomtype
                    else:
                        i_b = self.max_needed_atomtype
                else:
                    i_b = Intify(atomtype_i_tokens[1])

            else:
                i_a = i_b = Intify(atomtype_i_str)

            i_a_final = None
            i_b_final = None
            for i in range(i_a, i_b + 1):
                if ((i in self.needed_atomtypes) or (self.min_sel_atomtype <= i)):
                    i_a_final = i
                    break
            for i in reversed(range(i_a, i_b + 1)):
                if ((i in self.needed_atomtypes) or (self.max_sel_atomtype >= i)):
                    i_b_final = i
                    break

            if not (i_a_final and i_b_final and j_a_final and j_b_final):
                del self.l_in_masses[i_line]
            elif ('*' in atomtype_i_str):
                del self.l_in_masses[i_line]
                for i in range(i_a_final, i_b_final + 1):
                    tokens[1] = Stringify(i,
                                          self.atomtypes_int2name,
                                          '@','atom','type')
                    self.l_in_masses.insert(i_line, (' ' * self.indent) +
                                            (' '.join(tokens)))
                    i_line += 1
            else:
                assert(i_a == i_b)
                tokens[1] = Stringify(i_a,
                                      self.atomtypes_int2name,
                                      '@','atom','type')
                self.l_in_masses[i_line] = (' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # --- BONDS AND BOND COEFFS ---

        # delete lines from l_data_bonds if they involve atoms we don't care about
        i_line = 0
        while i_line < len(self.l_data_bonds):
            line = self.l_data_bonds[i_line]
            tokens = line.strip().split()
            assert(len(tokens) == 4)

            bondid = Intify(tokens[0])
            bondtype = Intify(tokens[1])
            atomid1 = Intify(tokens[2])
            atomid2 = Intify(tokens[3])
            # if ((atomid1 in self.needed_atomids) and
            #     (atomid2 in self.needed_atomids)):
            tokens[0] = '$bond:id' + str(bondid)
            tokens[1] = Stringify(bondtype,
                                  self.bondtypes_int2name,
                                  '@','bond','type')
            tokens[2] = Stringify(atomid1,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[3] = Stringify(atomid2,
                                  self.atomids_int2name,
                                  '$','atom','id')
            if self.ignore_bond_types:
                # Then instead of a "Bonds" section we want a "Bond List" 
                # section which omits the bond type information (in tokens[1])
                del tokens[1]
            self.needed_bondids.add(bondid)
            self.needed_bondtypes.add(bondtype)
            self.l_data_bonds[i_line] = (' ' * self.indent) + (' '.join(tokens))
            i_line += 1
            # else:
            #    del self.l_data_bonds[i_line]

        # delete data_bond_coeffs for bondtypes we no longer care about
        i_line = 0
        while i_line < len(self.l_data_bond_coeffs):
            line = self.l_data_bond_coeffs[i_line]
            tokens = line.strip().split()
            bondtype = Intify(tokens[0])
            if (not (bondtype in self.needed_bondtypes)):
                del self.l_data_bond_coeffs[i_line]
            else:
                tokens[0] = Stringify(bondtype,
                                      self.bondtypes_int2name,
                                      '@','bond','type')
                self.l_data_bond_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # delete in_bond_coeffs for bondtypes we no longer care about:
        if len(self.needed_bondtypes) > 0:
            for bondtype in self.needed_bondtypes:
                assert(type(bondtype) is int)
                if ((self.min_needed_bondtype == None) or
                    (self.min_needed_bondtype > bondtype)):
                    self.min_needed_bondtype = bondtype
                if ((self.max_needed_bondtype == None) or
                    (self.max_needed_bondtype < bondtype)):
                    self.max_needed_bondtype = bondtype
            for bondid in self.needed_bondids:
                assert(type(bondid) is int)
                if ((self.min_needed_bondid == None) or
                    (self.min_needed_bondid > bondid)):
                    self.min_needed_bondid = bondid
                if ((self.max_needed_bondid == None) or
                    (self.max_needed_bondid < bondid)):
                    self.max_needed_bondid = bondid
        else:
            # If no bonds were needed, then define some defaults
            # to make sure we don't keep any of them later.
            self.min_needed_bondtype = 1
            self.max_needed_bondtype = 0


        i_line = 0
        while i_line < len(self.l_in_bond_coeffs):
            line = self.l_in_bond_coeffs[i_line]
            tokens = line.strip().split()
            bondtype_str = tokens[1]

            if ('*' in bondtype_str):
                bondtype_tokens = bondtype_str.split('*')

                if bondtype_tokens[0] == '':
                    i_a = self.min_needed_bondtype
                else:
                    i_a = Intify(bondtype_tokens[0])

                if bondtype_tokens[1] == '':
                    i_b = self.max_needed_bondtype
                else:
                    i_b = Intify(bondtype_tokens[1])

            else:
                i_a = Intify(bondtype_str)
                i_b = i_a

            if i_a < self.min_needed_bondtype:
                i_a = self.min_needed_bondtype
            if i_b > self.max_needed_bondtype:
                i_b = self.max_needed_bondtype

            if ('*' in bondtype_str):
                del self.l_in_bond_coeffs[i_line]
                for i in range(i_a, i_b + 1):
                    if (i in self.needed_bondtypes):
                        tokens[1] = Stringify(i,
                                              self.bondtypes_int2name,
                                              '@','bond','type')
                        self.l_in_bond_coeffs.insert(i_line,
                                                     (' ' * self.indent) + (' '.join(tokens)))
                        i_line += 1
            else:
                if i_a < i_b:
                    raise InputError('Error: number of bond types in data file is not consistent with the\n'
                                     '       number of bond types you have define bond_coeffs for.\n')
                if (i_a == i_b) and (i_a in self.needed_bondtypes):
                    tokens[1] = Stringify(i_a,
                                          self.bondtypes_int2name,
                                          '@','bond','type')
                    self.l_in_bond_coeffs[i_line] = (
                        ' ' * self.indent) + (' '.join(tokens))
                    i_line += 1
                else:
                    del self.l_in_bond_coeffs[i_line]

        # --- ANGLES AND ANGLE COEFFS ---

        # delete lines from data_angles if they involve atoms we don't care about
        i_line = 0
        while i_line < len(self.l_data_angles):
            line = self.l_data_angles[i_line]
            tokens = line.strip().split()
            assert(len(tokens) == 5)

            angleid = Intify(tokens[0])
            angletype = Intify(tokens[1])
            atomid1 = Intify(tokens[2])
            atomid2 = Intify(tokens[3])
            atomid3 = Intify(tokens[4])
            # if ((atomid1 in self.needed_atomids) and
            #     (atomid2 in self.needed_atomids) and
            #     (atomid3 in self.needed_atomids)):
            tokens[0] = '$angle:id' + str(angleid)

            tokens[1] = Stringify(angletype,
                                  self.angletypes_int2name,
                                  '@','angle','type')
            tokens[2] = Stringify(atomid1,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[3] = Stringify(atomid2,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[4] = Stringify(atomid3,
                                  self.atomids_int2name,
                                  '$','atom','id')
            self.needed_angleids.add(angleid)
            self.needed_angletypes.add(angletype)
            self.l_data_angles[i_line] = (' ' * self.indent) + (' '.join(tokens))
            i_line += 1
            # else:
            #    del self.l_data_angles[i_line]

        # delete data_angle_coeffs for angletypes we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_angle_coeffs):
            line = self.l_data_angle_coeffs[i_line]
            tokens = line.strip().split()
            angletype = Intify(tokens[0])
            if (not (angletype in self.needed_angletypes)):
                del self.l_data_angle_coeffs[i_line]
            else:
                tokens[0] = Stringify(angletype,
                                      self.angletypes_int2name,
                                      '@','angle','type')
                self.l_data_angle_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # --- class2specific ----
        # Do the same for BondBond and BondAngle Coeffs:
        # NOTE: LAMMPS INPUT SCRIPTS, ALL CLASS2 COEFFS are represented by:
        #       angle_coeff, dihedral_coeff, and improper_coeff commands.
        #       THERE ARE NO bondbond_coeff commands, or bondangle_coeff commands,
        #       etc..., so we dont have to worry about l_in_bondbond_coeffs,...
        # Delete data_bondbond_coeffs for angletypes we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_bondbond_coeffs):
            line = self.l_data_bondbond_coeffs[i_line]
            tokens = line.strip().split()
            angletype = Intify(tokens[0])
            if (not (angletype in self.needed_angletypes)):
                del self.l_data_bondbond_coeffs[i_line]
            else:
                tokens[0] = Stringify(angletype,
                                      self.angletypes_int2name,
                                      '@','angle','type')
                self.l_data_bondbond_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # Delete data_bondangle_coeffs for angletypes we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_bondangle_coeffs):
            line = self.l_data_bondangle_coeffs[i_line]
            tokens = line.strip().split()
            angletype = Intify(tokens[0])
            if (not (angletype in self.needed_angletypes)):
                del self.l_data_bondangle_coeffs[i_line]
            else:
                tokens[0] = Stringify(angletype,
                                      self.angletypes_int2name,
                                      '@','angle','type')
                self.l_data_bondangle_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # --- end of class2specific ----

        # Delete in_angle_coeffs for angletypes we no longer care about:
        if len(self.needed_angletypes) > 0:
            for angletype in self.needed_angletypes:
                assert(type(angletype) is int)
                if ((self.min_needed_angletype == None) or
                        (self.min_needed_angletype > angletype)):
                    self.min_needed_angletype = angletype
                if ((self.max_needed_angletype == None) or
                        (self.max_needed_angletype < angletype)):
                    self.max_needed_angletype = angletype
            for angleid in self.needed_angleids:
                assert(type(angleid) is int)
                if ((self.min_needed_angleid == None) or
                        (self.min_needed_angleid > angleid)):
                    self.min_needed_angleid = angleid
                if ((self.max_needed_angleid == None) or
                        (self.max_needed_angleid < angleid)):
                    self.max_needed_angleid = angleid
        else:
            # If no angles were needed, then define some defaults
            # to make sure we don't keep any of them later.
            self.min_needed_angletype = 1
            self.max_needed_angletype = 0


        i_line = 0
        while i_line < len(self.l_in_angle_coeffs):
            line = self.l_in_angle_coeffs[i_line]
            tokens = line.strip().split()
            angletype_str = tokens[1]

            if ('*' in angletype_str):
                angletype_tokens = angletype_str.split('*')

                if angletype_tokens[0] == '':
                    i_a = self.min_needed_angletype
                else:
                    i_a = Intify(angletype_tokens[0])

                if angletype_tokens[1] == '':
                    i_b = self.max_needed_angletype
                else:
                    i_b = Intify(angletype_tokens[1])

            else:
                i_a = i_b = Intify(angletype_str)

            if i_a < self.min_needed_angletype:
                i_a = self.min_needed_angletype
            if i_b > self.max_needed_angletype:
                i_b = self.max_needed_angletype

            if ('*' in angletype_str):
                del self.l_in_angle_coeffs[i_line]
                for i in range(i_a, i_b + 1):
                    if (i in self.needed_angletypes):
                        tokens[1] = Stringify(i,
                                              self.angletypes_int2name,
                                              '@','angle','type')
                        self.l_in_angle_coeffs.insert(i_line,
                                                      (' ' * self.indent) + (' '.join(tokens)))
                        i_line += 1
            else:
                if i_a < i_b:
                    raise InputError('Error: number of angle types in data file is not consistent with the\n'
                                     '       number of angle types you have define angle_coeffs for.\n')
                if (i_a == i_b) and (i_a in self.needed_angletypes):
                    tokens[1] = Stringify(i_a,
                                          self.angletypes_int2name,
                                          '@','angle','type')
                    self.l_in_angle_coeffs[i_line] = (
                        ' ' * self.indent) + (' '.join(tokens))
                    i_line += 1
                else:
                    del self.l_in_angle_coeffs[i_line]

        # --- DIHEDRALS AND DIHEDRAL COEFFS ---

        # delete lines from data_dihedrals if they involve atoms we don't care
        # about
        i_line = 0
        while i_line < len(self.l_data_dihedrals):
            line = self.l_data_dihedrals[i_line]
            tokens = line.strip().split()
            assert(len(tokens) == 6)

            dihedralid = Intify(tokens[0])
            dihedraltype = Intify(tokens[1])
            atomid1 = Intify(tokens[2])
            atomid2 = Intify(tokens[3])
            atomid3 = Intify(tokens[4])
            atomid4 = Intify(tokens[5])

            # if ((atomid1 in self.needed_atomids) and
            #     (atomid2 in self.needed_atomids) and
            #     (atomid3 in self.needed_atomids) and
            #     (atomid4 in self.needed_atomids)):

            tokens[0] = '$dihedral:id' + str(dihedralid)
            tokens[1] = Stringify(dihedraltype,
                                  self.dihtypes_int2name,
                                  '@','dihedral','type')
            tokens[2] = Stringify(atomid1,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[3] = Stringify(atomid2,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[4] = Stringify(atomid3,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[5] = Stringify(atomid4,
                                  self.atomids_int2name,
                                  '$','atom','id')
            self.needed_dihedralids.add(dihedralid)
            self.needed_dihedraltypes.add(dihedraltype)
            self.l_data_dihedrals[i_line] = (' ' * self.indent) + (' '.join(tokens))
            i_line += 1
            # else:
            #    del self.l_data_dihedrals[i_line]

        # delete data_dihedral_coeffs for dihedraltypes we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_dihedral_coeffs):
            line = self.l_data_dihedral_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_dihedral_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_dihedral_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # --- class2specific ----
        # Do the same for MiddleBondTorsion, EndBondTorsion, AngleTorsion,
        #                 AngleAngleTorsion, and BondBond13 Coeffs
        # NOTE: LAMMPS INPUT SCRIPTS, ALL CLASS2 COEFFS are represented by:
        #       angle_coeff, dihedral_coeff, and improper_coeff commands.
        #       THERE ARE NO "middlebondtorsion_coeff" commands, etc...so we don't
        #       have to worry about dealing with "self.l_in_middlebondtorsion_coeffs",...
        # delete data_middlebondtorsion_coeffs for dihedraltypes
        # we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_middlebondtorsion_coeffs):
            line = self.l_data_middlebondtorsion_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_middlebondtorsion_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_middlebondtorsion_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # delete data_endbondtorsion_coeffs for dihedraltypes we 
        # no longer care about:
        i_line = 0
        while i_line < len(self.l_data_endbondtorsion_coeffs):
            line = self.l_data_endbondtorsion_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_endbondtorsion_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_endbondtorsion_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # delete data_angletorsion_coeffs for dihedraltypes we 
        # no longer care about:
        i_line = 0
        while i_line < len(self.l_data_angletorsion_coeffs):
            line = self.l_data_angletorsion_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_angletorsion_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_angletorsion_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # delete data_angleangletorsion_coeffs for dihedraltypes we 
        # no longer care about:
        i_line = 0
        while i_line < len(self.l_data_angleangletorsion_coeffs):
            line = self.l_data_angleangletorsion_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_angleangletorsion_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_angleangletorsion_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # delete data_bondbond13_coeffs for dihedraltypes we 
        # no longer care about:
        i_line = 0
        while i_line < len(self.l_data_bondbond13_coeffs):
            line = self.l_data_bondbond13_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype = Intify(tokens[0])
            if (not (dihedraltype in self.needed_dihedraltypes)):
                del self.l_data_bondbond13_coeffs[i_line]
            else:
                tokens[0] = Stringify(dihedraltype,
                                      self.dihtypes_int2name,
                                      '@','dihedral','type')
                self.l_data_bondbond13_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # --- end of class2specific ----


        # delete in_dihedral_coeffs for dihedraltypes we no longer care about:
        if len(self.needed_dihedraltypes) > 0:
            for dihedraltype in self.needed_dihedraltypes:
                assert(type(dihedraltype) is int)
                if ((self.min_needed_dihedraltype == None) or
                        (self.min_needed_dihedraltype > dihedraltype)):
                    self.min_needed_dihedraltype = dihedraltype
                if ((self.max_needed_dihedraltype == None) or
                        (self.max_needed_dihedraltype < dihedraltype)):
                    self.max_needed_dihedraltype = dihedraltype
            for dihedralid in self.needed_dihedralids:
                assert(type(dihedralid) is int)
                if ((self.min_needed_dihedralid == None) or
                        (self.min_needed_dihedralid > dihedralid)):
                    self.min_needed_dihedralid = dihedralid
                if ((self.max_needed_dihedralid == None) or
                        (self.max_needed_dihedralid < dihedralid)):
                    self.max_needed_dihedralid = dihedralid
        else:
            # If no dihedrals were needed, then define some defaults
            # to make sure we don't keep any of them later.
            self.min_needed_dihedraltype = 1
            self.max_needed_dihedraltype = 0


        i_line = 0
        while i_line < len(self.l_in_dihedral_coeffs):
            line = self.l_in_dihedral_coeffs[i_line]
            tokens = line.strip().split()
            dihedraltype_str = tokens[1]

            if ('*' in dihedraltype_str):
                dihedraltype_tokens = dihedraltype_str.split('*')

                if dihedraltype_tokens[0] == '':
                    i_a = self.min_needed_dihedraltype
                else:
                    i_a = Intify(dihedraltype_tokens[0])

                if dihedraltype_tokens[1] == '':
                    i_b = self.max_needed_dihedraltype
                else:
                    i_b = Intify(dihedraltype_tokens[1])

            else:
                i_a = i_b = Intify(dihedraltype_str)

            if i_a < self.min_needed_dihedraltype:
                i_a = self.min_needed_dihedraltype
            if i_b > self.max_needed_dihedraltype:
                i_b = self.max_needed_dihedraltype

            if ('*' in dihedraltype_str):
                del self.l_in_dihedral_coeffs[i_line]
                for i in range(i_a, i_b + 1):
                    if (i in self.needed_dihedraltypes):
                        tokens[1] = Stringify(i,
                                              self.dihtypes_int2name,
                                              '@','dihedral','type')
                        self.l_in_dihedral_coeffs.insert(i_line,
                                                         (' ' * self.indent) + (' '.join(tokens)))
                        i_line += 1
            else:
                if i_a < i_b:
                    raise InputError('Error: number of dihedral types in data file is not consistent with the\n'
                                     '       number of dihedral types you have define dihedral_coeffs for.\n')
                if (i_a == i_b) and (i_a in self.needed_dihedraltypes):
                    tokens[1] = Stringify(i_a,
                                          self.dihtypes_int2name,
                                          '@','dihedral','type')
                    self.l_in_dihedral_coeffs[i_line] = (
                        ' ' * self.indent) + (' '.join(tokens))
                    i_line += 1
                else:
                    del self.l_in_dihedral_coeffs[i_line]

        # --- IMPROPERS AND IMPROPER COEFFS ---

        # delete lines from data_impropers if they involve atoms we don't care
        # about
        i_line = 0
        while i_line < len(self.l_data_impropers):
            line = self.l_data_impropers[i_line]
            tokens = line.strip().split()
            assert(len(tokens) == 6)

            improperid = Intify(tokens[0])
            impropertype = Intify(tokens[1])
            atomid1 = Intify(tokens[2])
            atomid2 = Intify(tokens[3])
            atomid3 = Intify(tokens[4])
            atomid4 = Intify(tokens[5])

            tokens[0] = '$improper:id' + str(improperid)
            tokens[1] = '@improper:type' + str(impropertype)

            tokens[1] = Stringify(impropertype,
                                  self.imptypes_int2name,
                                  '@','improper','type')
            tokens[2] = Stringify(atomid1,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[3] = Stringify(atomid2,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[4] = Stringify(atomid3,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[5] = Stringify(atomid4,
                                  self.atomids_int2name,
                                  '$','atom','id')
            self.needed_improperids.add(improperid)
            self.needed_impropertypes.add(impropertype)
            self.l_data_impropers[i_line] = (' ' * self.indent) + (' '.join(tokens))
            i_line += 1
            # else:
            #    del self.l_data_impropers[i_line]


        # delete data_improper_coeffs for impropertypes we no longer care about:
        i_line = 0
        while i_line < len(self.l_data_improper_coeffs):
            line = self.l_data_improper_coeffs[i_line]
            tokens = line.strip().split()
            impropertype = Intify(tokens[0])
            if (not (impropertype in self.needed_impropertypes)):
                del self.l_data_improper_coeffs[i_line]
            else:
                tokens[0] = Stringify(impropertype,
                                      self.imptypes_int2name,
                                      '@','improper','type')
                self.l_data_improper_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1

        # --- class2specific ----
        # Do the same for AngleAngle Coeffs
        # NOTE: LAMMPS INPUT SCRIPTS, ALL CLASS2 COEFFS are represented by:
        #       angle_coeff, dihedral_coeff, and improper_coeff commands.
        #       THERE ARE NO "angleangle_coeff" commands, etc...so we don't
        #       have to worry about dealing with "l_in_angleangle_coeffs",...
        # delete entries in l_data_angleangle_coeffs for impropertypes we
        # no longer care about:
        i_line = 0
        while i_line < len(self.l_data_angleangle_coeffs):
            line = self.l_data_angleangle_coeffs[i_line]
            tokens = line.strip().split()
            impropertype = Intify(tokens[0])
            if (not (impropertype in self.needed_impropertypes)):
                del self.l_data_angleangle_coeffs[i_line]
            else:
                tokens[0] = Stringify(impropertype,
                                      self.imptypes_int2name,
                                      '@','improper','type')
                self.l_data_angleangle_coeffs[i_line] = (
                    ' ' * self.indent) + (' '.join(tokens))
                i_line += 1
        # --- end of class2specific ----

        # delete in_improper_coeffs for impropertypes we no longer care about:
        if len(self.needed_impropertypes) > 0:
            for impropertype in self.needed_impropertypes:
                assert(type(impropertype) is int)
                if ((self.min_needed_impropertype == None) or
                    (self.min_needed_impropertype > impropertype)):
                    self.min_needed_impropertype = impropertype
                if ((self.max_needed_impropertype == None) or
                    (self.max_needed_impropertype < impropertype)):
                    self.max_needed_impropertype = impropertype
            for improperid in self.needed_improperids:
                assert(type(improperid) is int)
                if ((self.min_needed_improperid == None) or
                    (self.min_needed_improperid > improperid)):
                    self.min_needed_improperid = improperid
                if ((self.max_needed_improperid == None) or
                    (self.max_needed_improperid < improperid)):
                    self.max_needed_improperid = improperid
        else:
            # If no impropers were needed, then define some defaults
            # to make sure we don't keep any of them later.
            self.min_needed_impropertype = 1
            self.max_needed_impropertype = 0


        i_line = 0
        while i_line < len(self.l_in_improper_coeffs):
            line = self.l_in_improper_coeffs[i_line]
            tokens = line.strip().split()
            impropertype_str = tokens[1]

            if ('*' in impropertype_str):
                impropertype_tokens = impropertype_str.split('*')

                if impropertype_tokens[0] == '':
                    i_a = self.min_needed_impropertype
                else:
                    i_a = Intify(impropertype_tokens[0])

                if impropertype_tokens[1] == '':
                    i_b = self.max_needed_impropertype
                else:
                    i_b = Intify(impropertype_tokens[1])

            else:
                i_a = i_b = Intify(impropertype_str)

            assert((type(i_a) is int) and (type(i_b) is int))

            if (i_a < self.min_needed_impropertype):
                i_a = self.min_needed_impropertype
            if (i_b > self.max_needed_impropertype):
                i_b = self.max_needed_impropertype

            if ('*' in impropertype_str):
                del self.l_in_improper_coeffs[i_line]
                for i in range(i_a, i_b + 1):
                    if (i in self.needed_impropertypes):
                        tokens[1] = Stringify(i,
                                              self.imptypes_int2name,
                                              '@','improper','type')
                        self.l_in_improper_coeffs.insert(i_line,
                                                         (' ' * self.indent) + (' '.join(tokens)))
                        i_line += 1
            else:
                if i_a < i_b:
                    raise InputError('Error: number of improper types in data file is not consistent with the\n'
                                     '       number of improper types you have define improper_coeffs for.\n')
                if (i_a == i_b) and (i_a in self.needed_impropertypes):
                    tokens[1] = Stringify(i_a,
                                          self.imptypes_int2name,
                                          '@','improper','type')
                    self.l_in_improper_coeffs[i_line] = (
                        ' ' * self.indent) + (' '.join(tokens))
                    i_line += 1
                else:
                    del self.l_in_improper_coeffs[i_line]



        # --- CMAP INTERACTIONS ---

        # Add the correct variable types to each line in self.l_data_cmaps
        # ($cmap: $atom:)
        # Also: delete lines from data_cmap if they involve atoms we don't care
        # about
        i_line = 0
        while i_line < len(self.l_data_cmap):
            line = self.l_data_cmap[i_line]
            tokens = line.strip().split()
            assert(len(tokens) == 7)

            cmapid = Intify(tokens[0])
            cmaptype = Intify(tokens[1])
            atomid1 = Intify(tokens[2])
            atomid2 = Intify(tokens[3])
            atomid3 = Intify(tokens[4])
            atomid4 = Intify(tokens[5])
            atomid5 = Intify(tokens[6])
            # if ((atomid1 in self.needed_atomids) and
            #    (atomid2 in self.needed_atomids)):
            tokens[0] = '$cmap:id' + str(cmapid)
            #tokens[1] = '@cmap:type' + str(cmaptype)
            tokens[1] = str(cmaptype)
            tokens[2] = Stringify(atomid1,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[3] = Stringify(atomid2,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[4] = Stringify(atomid3,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[5] = Stringify(atomid4,
                                  self.atomids_int2name,
                                  '$','atom','id')
            tokens[5] = Stringify(atomid5,
                                  self.atomids_int2name,
                                  '$','atom','id')


            self.needed_cmapids.add(cmapid)
            self.needed_cmaptypes.add(cmaptype)
            self.l_data_cmap[i_line] = (' ' * self.indent) + (' '.join(tokens))
            i_line += 1
            # else:
            #    del self.l_data_cmap[i_line]

        ## delete in_cmap_coeffs for cmaptypes we no longer care about:
        #if len(self.needed_cmaptypes) > 0:
        #    for cmaptype in self.needed_cmaptypes:
        #        assert(type(cmaptype) is int)
        #        if ((self.min_needed_cmaptype == None) or
        #            (self.min_needed_cmaptype > cmaptype)):
        #            self.min_needed_cmaptype = cmaptype
        #        if ((self.max_needed_cmaptype == None) or
        #            (self.max_needed_cmaptype < cmaptype)):
        #            self.max_needed_cmaptype = cmaptype
        #    for cmapid in self.needed_cmapids:
        #        assert(type(cmapid) is int)
        #        if ((self.min_needed_cmapid == None) or
        #            (self.min_needed_cmapid > cmapid)):
        #            self.min_needed_cmapid = cmapid
        #        if ((self.max_needed_cmapid == None) or
        #            (self.max_needed_cmapid < cmapid)):
        #            self.max_needed_cmapid = cmapid
        #else:
        #    # If no cmap interactions were needed, then define some defaults
        #    # to make sure we don't keep any of them later.
        #    self.min_needed_cmaptype = 1
        #    self.max_needed_cmaptype = 0
        #
        #
        #i_line = 0
        #while i_line < len(self.l_in_cmap_coeffs):
        #    line = self.l_in_cmap_coeffs[i_line]
        #    tokens = line.strip().split()
        #    cmaptype_str = tokens[1]
        #
        #    if ('*' in cmaptype_str):
        #        cmaptype_tokens = cmaptype_str.split('*')
        #
        #        if cmaptype_tokens[0] == '':
        #            i_a = self.min_needed_cmaptype
        #        else:
        #            i_a = Intify(cmaptype_tokens[0])
        #
        #        if cmaptype_tokens[1] == '':
        #            i_b = self.max_needed_cmaptype
        #        else:
        #            i_b = Intify(cmaptype_tokens[1])
        #
        #    else:
        #        i_a = i_b = Intify(cmaptype_str)
        #
        #    assert((type(i_a) is int) and (type(i_b) is int))
        #
        #    if (i_a < self.min_needed_cmaptype):
        #        i_a = self.min_needed_cmaptype
        #    if (i_b > self.max_needed_cmaptype):
        #        i_b = self.max_needed_cmaptype
        #    
        #    # if i_a == i_b:
        #         tokens[1] = Stringify(atomid5,
        #                               self.cmaptypes_int2name,
        #                               '@','cmap','type')
        #    
        #    if ('*' in cmaptype_str):
        #        del self.l_in_cmap_coeffs[i_line]
        #        for i in range(i_a, i_b + 1):
        #            if (i in self.needed_cmaptypes):
        #               tokens[1] = Stringify(i,
        #                                     self.cmaptypes_int2name,
        #                                     '@','cmap','type')
        #                self.l_in_cmap_coeffs.insert(i_line,
        #                                             (' ' * self.indent) + (' '.join(tokens)))
        #                i_line += 1
        #    else:
        #        if i_a < i_b:
        #            raise InputError('Error: number of cmap types in data file is not consistent with the\n'
        #                             '       number of cmap types you have define cmap_coeffs for.\n')
        #        if (i_a == i_b) and (i_a in self.needed_cmaptypes):
        #            tokens[1] = Stringify(i_a,
        #                                  self.cmaptypes_int2name,
        #                                  '@','cmap','type')
        #            self.l_in_cmap_coeffs[i_line] = (
        #                ' ' * self.indent) + (' '.join(tokens))
        #            i_line += 1
        #        else:
        #            del self.l_in_cmap_coeffs[i_line]





        # --- GROUPS ---

        # Now parse through all of the "group" commands and try and figure
        # out if any of these groups contain any of the atoms we are keeping.
        # If so, then save the group and write it out.
        # (I hate trying to parse this kind of text.)

        # if len(self.l_in_group) > 0:
        #    sys.stderr.write('\n'
        #                     ' --groups--  Attempting to parse \"group\" commands.\n'
        #                     '         This may cause '+g_program_name+' to crash.\n'
        #                     '         If so, comment out all group commands in your input script(s), and\n'
        #                     '         try again.  (And please report the error. -Andrew 2017-10)\n')

        i_line = 0
        groups_needed = set(['all'])
        while i_line < len(self.l_in_group):
            line = self.l_in_group[i_line]
            tokens = line.strip().split()
            delete_this_command = False
            explicit_definition = False
            if len(tokens) < 3:
                delete_this_command = True
            group_name = tokens[1]
            specifier_style = tokens[2]
            str_logical = ''
            str_selection = ''
            if specifier_style[0:4] == 'type':
                str_logical += specifier_style[4:]
                explicit_definition = True
                specifier_style = 'type'
            elif specifier_style == 'id':
                str_logical += specifier_style[2:]
                explicit_definition = True
                specifier_style = 'id'
            elif specifier_style == 'molecule':
                str_logical += specifier_style[8:]
                specifier_style = 'molecule'
                explicit_definition = True

            if explicit_definition:
                i_token_sel_min = 3
                if len(tokens) <= i_token_sel_min:
                    sys.stderr.write('WARNING: possible syntax error on this line:\n'
                                     + '        ' + self.l_in_group[i_line] + '\n')
                    delete_this_command = True
                if str_logical == '':
                    str_logical = tokens[i_token_sel_min]
                    if not str_logical[0].isdigit():
                        i_token_sel_min += 1
                        if len(tokens) <= i_token_sel_min:
                            tokens.append('')
                else:
                    tokens.insert(i_token_sel_min, str_logical)

                i_token_sel_max = len(tokens) - 1

                for i in range(i_token_sel_min, len(tokens)):
                    if tokens[i].isdigit():
                        break
                    else:
                        i_token_sel_max = i

                assert(len(tokens) > i_token_sel_min)

                if str_logical[0:2] in ('<=', '>=', '==', '!=', '<>'):
                    tokens[i_token_sel_min] = str_logical[
                        2:] + tokens[i_token_sel_min]
                    str_logical = str_logical[0:2]
                    if str_logical == '<=':
                        self.l_group_selection = [(None, int(tokens[i_token_sel_min]))]
                    elif str_logical == '>=':
                        self.l_group_selection = [(int(tokens[i_token_sel_min]), None)]
                    elif str_logical == '==':
                        self.l_group_selection = [(int(tokens[i_token_sel_min]),
                                              int(tokens[i_token_sel_min]))]
                    elif str_logical == '!=':
                        self.l_group_selection = [(None, int(tokens[i_token_sel_min]) - 1),
                                             (int(tokens[i_token_sel_min]) + 1, None)]
                    elif str_logical == '<>':
                        self.l_group_selection= [(int(tokens[i_token_sel_min]),
                                                  int(tokens[i_token_sel_max]))]

                elif str_logical[0:1] in ('<', '>'):
                    tokens[i_token_sel_min] = str_logical[
                        1:] + tokens[i_token_sel_min]
                    str_logical = str_logical[0:1]
                    if str_logical == '<':
                        self.l_group_selection = [
                            (None, int(tokens[i_token_sel_min]) - 1)]
                    elif str_logical == '>':
                        self.l_group_selection = [
                            (int(tokens[i_token_sel_min]) + 1, None)]
                else:
                    str_selection = ' '.join(
                        tokens[i_token_sel_min:i_token_sel_max + 1])
                    self.l_group_selection = LammpsSelectToIntervals(str_selection,
                                                                     slice_delim=':',
                                                                     or_delim=' ')

                mn, mx = IntervalListToMinMax(self.l_group_selection)
                if mn == None:
                    mn = 1
                filtered_selection = []
                if specifier_style == 'type':
                    if mx == None:
                        mx = self.max_needed_atomtype
                    for i in range(mn, mx + 1):
                        if (BelongsToSel(i, self.l_group_selection)
                                and (i in self.needed_atomtypes)):
                            filtered_selection.append((i, i))
                elif specifier_style == 'id':
                    if mx == None:
                        mx = self.max_needed_atomid
                    for i in range(mn, mx + 1):
                        if (BelongsToSel(i, self.l_group_selection)
                                and (i in self.needed_atomids)):
                            filtered_selection.append((i, i))
                elif specifier_style == 'molecule':
                    if mx == None:
                        mx = self.max_needed_molid
                    for i in range(mn, mx + 1):
                        if (BelongsToSel(i, self.l_group_selection)
                                and (i in self.needed_molids)):
                            filtered_selection.append((i, i))

                MergeIntervals(filtered_selection)

                if len(filtered_selection) > 0:

                    tokens = ['group', group_name, specifier_style]
                    for interval in filtered_selection:
                        a = interval[0]
                        b = interval[1]

                        if specifier_style == 'type':
                            if a == b:
                                var_descr = Stringify(a,
                                                      self.atomtypes_int2name,
                                                      '@','atom','type')
                                tokens.append(var_descr)
                            else:
                                var_name_a = LookupVarName(a,
                                                           self.atomtypes_int2name,
                                                           'type')
                                var_name_b = LookupVarName(b,
                                                           self.atomtypes_int2name,
                                                           'type')
                                tokens.append('@{atom:' + var_name_a +
                                              '}:@{atom:' + var_name_b + '}')

                        if specifier_style == 'id':
                            if a == b:
                                var_descr = Stringify(a,
                                                      self.atomids_int2name,
                                                      '$','atom','id')
                                tokens.append(var_descr)
                            else:
                                var_name_a = LookupVarName(a,
                                                           self.atomids_int2name,
                                                           'id')
                                var_name_b = LookupVarName(b,
                                                           self.atomids_int2name,
                                                           'id')
                                tokens.append('${atom:' + var_name_a +
                                              '}:${atom:' + var_name_b + '}')


                        if specifier_style == 'molecule':
                            if a == b:
                                var_descr = Stringify(a,
                                                      self.molids_int2name,
                                                      '$','mol','id')
                                tokens.append(var_descr)
                            else:
                                var_name_a = LookupVarName(a,
                                                           self.molids_int2name,
                                                           'id')
                                var_name_b = LookupVarName(b,
                                                           self.molids_int2name,
                                                           'id')
                                tokens.append('${atom:' + var_name_a +
                                              '}:${atom:' + var_name_b + '}')


                    # Commenting out next two lines.  (This is handled later.)
                    #self.l_in_group[i_line] = ' '.join(tokens)
                    # groups_needed.add(group_name)

                else:
                    delete_this_command = True

            else:
                if len(tokens) > 3:
                    if tokens[2] == 'union':
                        i_token = 3
                        while i_token < len(tokens):
                            if not (tokens[i_token] in groups_needed):
                                del tokens[i_token]
                            else:
                                i_token += 1
                        # if none of the groups contain atoms we need,
                        # then delete the entire command
                        if len(tokens) <= 3:
                            delete_this_command = True
                    elif tokens[2] == 'intersect':
                        i_token = 3
                        while i_token < len(tokens):
                            if not (tokens[i_token] in groups_needed):
                                # if any of the groups we need are empty
                                # then delete the command
                                delete_this_command = True
                                break
                            i_token += 1
                    elif (tokens[2] == 'subtract') and (len(tokens) >= 5):
                        if not (tokens[3] in groups_needed):
                            delete_this_command = True
                        i_token = 4
                        while i_token < len(tokens):
                            if not (tokens[i_token] in groups_needed):
                                del tokens[i_token]
                            else:
                                i_token += 1
                    else:
                        # Otherwise I don't recongize the syntax of this
                        # group command.  In that case, I just delete it.
                        delete_this_command = True

                elif tokens[2] == 'clear':
                    pass
                elif tokens[2] == 'delete':
                    pass
                else:
                    delete_this_command = True
            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_group[i_line].rstrip() + '\"\n')
                del self.l_in_group[i_line]
            else:
                groups_needed.add(group_name)
                self.l_in_group[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1

        # --- fix rigid ---

        i_line = 0
        while i_line < len(self.l_in_fix_rigid):
            line = self.l_in_fix_rigid[i_line]
            tokens = line.strip().split()
            if len(tokens) < 4:
                break
            fixid = tokens[1]
            group_name = tokens[2]
            delete_this_command = True
            assert(tokens[3].find('rigid') == 0)
            if group_name in groups_needed:
                delete_this_command = False

            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_fix_rigid[i_line].rstrip() + '\"\n')
                del self.l_in_fix_rigid[i_line]
            else:
                self.l_in_fix_rigid[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1

        # --- set ---

        i_line = 0
        while i_line < len(self.l_in_set):
            line = self.l_in_set[i_line]
            tokens = line.strip().split()
            self.l_new_set_commands = []
            self.l_new_set_static_commands = []
            if len(tokens) < 4:
                break
            if tokens[1] == 'type':
                pattern = tokens[2].split('*')
                if pattern[0] == '':
                    types_lo = self.min_needed_atomtype
                else:
                    types_lo = types_hi = int(pattern[0])
                    if types_lo < self.min_needed_atomtype:
                        types_lo = self.min_needed_atomtype
                if len(pattern)  == 2:
                    if pattern[1] == '':
                        types_hi = self.max_needed_atomtype
                    else:
                        types_hi = min(int(pattern[1]), self.max_needed_atomtype)
                for i in range(types_lo, types_hi+1):
                    if i in self.needed_atomtypes:
                        var_descr = Stringify(i,
                                              self.atomtypes_int2name,
                                              '@','atom','type')
                        self.l_new_set_static_commands.append((' ' * self.indent) +
                                                              ' '.join(tokens[0:2])+' '+
                                                              var_descr + ' ' +
                                                              ' '.join(tokens[3:]))
            elif tokens[1] == 'atom':
                pattern = tokens[2].split('*')
                if pattern[0] == '':
                    atomids_lo = self.min_needed_atomid
                else:
                    atomids_lo = atomids_hi = int(pattern[0])
                    if atomids_lo < self.min_needed_atomid:
                        atomids_lo = self.min_needed_atomid
                if len(pattern)  == 2:
                    if pattern[1] == '':
                        atomids_hi = self.max_needed_atomid
                    else:
                        atomids_hi = min(int(pattern[1]), self.max_needed_atomid)
                for i in range(atomids_lo, atomids_hi+1):
                    if i in self.needed_atomids:
                        var_descr = Stringify(i,
                                              self.atomids_int2name,
                                              '$','atom','id')
                        self.l_new_set_commands.append((' ' * self.indent) +
                                                       ' '.join(tokens[0:2])+' '+
                                                       var_descr + ' ' +
                                                       ' '.join(tokens[3:]))
            elif tokens[1] == 'mol':
                pattern = tokens[2].split('*')
                if pattern[0] == '':
                    molids_lo = self.min_needed_molid
                else:
                    molids_lo = molids_hi = int(pattern[0])
                    if molids_lo < self.min_needed_molid:
                        molids_lo = self.min_needed_molid
                if len(pattern)  == 2:
                    if pattern[1] == '':
                        molids_hi = self.max_needed_molid
                    else:
                        molids_hi = min(int(pattern[1]), self.max_needed_molid)
                for i in range(molids_lo, molids_hi+1):
                    if i in self.needed_molids:
                        var_descr = Stringify(i,
                                              self.molids_int2name,
                                              '$','mol','id')
                        self.l_new_set_commands.append(' '.join(tokens[0:2])+' '+
                                                       var_descr + ' ' +
                                                       ' '.join(tokens[3:]))
            elif tokens[0] == 'group':
                group_name = tokens[2]
                if group_name in groups_needed:
                    self.l_new_set_static_commands = [self.l_in_set[i_line]]

            if len(self.l_new_set_commands) > 0:
                self.l_in_set[i_line:i_line+1] = self.l_new_set_commands
                i_line += len(self.l_new_set_commands)
            elif len(self.l_new_set_static_commands) > 0:
                self.l_in_set_static += self.l_new_set_static_commands
                del self.l_in_set[i_line]
            else:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_set[i_line].rstrip() + '\"\n')
                del self.l_in_set[i_line]
            

        # --- fix shake ---

        i_line = 0
        while i_line < len(self.l_in_fix_shake):
            line = self.l_in_fix_shake[i_line]
            tokens = line.strip().split()
            if len(tokens) < 4:
                break
            fixid = tokens[1]
            group_name = tokens[2]
            delete_this_command = True
            assert(tokens[3].find('shake') == 0)

            #  parse the list of angle types
            #i_token = tokens.index('a')
            for i_token in range(0, len(tokens)):
                if tokens[i_token] == 'a':
                    break
            if i_token != len(tokens):
                i_token += 1
                while (i_token < len(tokens)) and tokens[i_token].isdigit():
                    # delete angle types from the list which
                    # do not belong to the selection
                    btype = int(tokens[i_token])
                    if int(tokens[i_token]) in self.needed_angletypes:
                        var_descr = Stringify(btype,
                                              self.angletypes_int2name,
                                              '@','angle','type')
                        tokens[i_token] = var_descr
                        i_token += 1
                        delete_this_command = False
                    else:
                        del tokens[i_token]

            #  parse the list of bond types
            #i_token = tokens.index('b')
            for i_token in range(0, len(tokens)):
                if tokens[i_token] == 'b':
                    break
            if i_token != len(tokens):
                i_token += 1
                while (i_token < len(tokens)) and tokens[i_token].isdigit():
                    # delete bond types from the list which
                    # do not belong to the selection
                    btype = int(tokens[i_token])
                    if int(tokens[i_token]) in self.needed_bondtypes:
                        var_descr = Stringify(btype,
                                              self.bondtypes_int2name,
                                              '@','bond','type')
                        tokens[i_token] = var_descr
                        i_token += 1
                        delete_this_command = False
                    else:
                        del tokens[i_token]

            #  parse the list of atom types
            # i_token = tokens.index('t')
            for i_token in range(0, len(tokens)):
                if tokens[i_token] == 't':
                    break
            if i_token != len(tokens):
                i_token += 1
                while (i_token < len(tokens)) and tokens[i_token].isdigit():
                    # delete atom types from the list which
                    # do not belong to the selection
                    btype = int(tokens[i_token])
                    if int(tokens[i_token]) in self.needed_atomtypes:
                        var_descr = Stringify(btype,
                                              self.atomtypes_int2name,
                                              '@','atom','type')
                        tokens[i_token] = var_descr
                        i_token += 1
                        delete_this_command = False
                    else:
                        del tokens[i_token]

            #  Selecting atoms by mass feature should still work, so we
            #  don't need to delete or ignore these kinds of commands.
            # for i_token in range(0, len(tokens)):
            #    if tokens[i_token] == 'm':
            #        break
            # if i_token != len(tokens):
            #    delete_this_command = True

            if 'mol' in tokens:
                # What does the 'mol' keyword do?
                #               https://lammps.sandia.gov/doc/fix_shake.html
                #               (Excerpt below:)
                #
                # mol value = template-ID
                # template-ID = ID of molecule template specified
                #               in a separate molecule command.  See:
                #
                # ltemplify.py does not yet know how to parse files read by the
                # molecule command (https://lammps.sandia.gov/doc/molecule.html)
                # so I ignore these commands for now.  -Andrew 2019-11-11
                delete_this_command = True

            if not (group_name in groups_needed):
                delete_this_command = True

            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_fix_shake[i_line].rstrip() + '\"\n')
                del self.l_in_fix_shake[i_line]
            else:
                self.l_in_fix_shake[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1

        # --- fix poems ---

        i_line = 0
        while i_line < len(self.l_in_fix_poems):
            line = self.l_in_fix_poems[i_line]
            tokens = line.strip().split()
            if len(tokens) < 4:
                break
            fixid = tokens[1]
            group_name = tokens[2]
            delete_this_command = True
            assert(tokens[3].find('poems') == 0)
            if group_name in groups_needed:
                delete_this_command = False
            if tokens[4] != 'molecule':
                delete_this_command = True
                sys.stderr.write('WARNING: ' + g_program_name + ' ONLY supports \"fix poems\" commands\n'
                                 '         which use the \"molecule\" keyword.\n')
            if tokens[4] == 'file':
                sys.stderr.write('         If you want use external files with fix poems, then you will have to\n'
                                 '         generate the file yourself.  You ask use moltemplate to generate\n'
                                 '         this file for you, by manually adding a section at the end of your\n'
                                 '         final .LT file (eg. \"system.lt\") which resembles the following:\n\n'
                                 'write(\"poems_file.txt\") {\n'
                                 '  1 1 $atom:idname1a $atom:idname2a $atom:idname3a ...\n'
                                 '  2 1 $atom:idname1b $atom:idname2b $atom:idname3b ...\n'
                                 '  3 1 $atom:idname1c $atom:idname2c $atom:idname3c ...\n'
                                 '  : :   etc...\n'
                                 '}\n\n'
                                 '      ...where $atom:idname1a, $atom:idname2a, ... are moltemplate-compatible\n'
                                 '         unique (full,long) id-names for the atoms in each rigid body.\n'
                                 '         This will insure the atom-id numbers in this file are correct.\n'

                                 '         See the documentation for fix poems for details.\n')

            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_fix_poems[i_line].rstrip() + '\"\n')
                del self.l_in_fix_poems[i_line]
            else:
                self.l_in_fix_poems[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1

        # --- fix qeq ---

        i_line = 0
        while i_line < len(self.l_in_fix_qeq):
            line = self.l_in_fix_qeq[i_line]
            tokens = line.strip().split()
            if len(tokens) < 4:
                break
            fixid = tokens[1]
            group_name = tokens[2]
            delete_this_command = True
            assert(tokens[3].find('qeq') == 0)
            if group_name in groups_needed:
                delete_this_command = False

            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_fix_qeq[i_line].rstrip() + '\"\n')
                del self.l_in_fix_qeq[i_line]
            else:
                self.l_in_fix_qeq[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1

        # --- fix qmmm ---

        i_line = 0
        while i_line < len(self.l_in_fix_qmmm):
            line = self.l_in_fix_qmmm[i_line]
            tokens = line.strip().split()
            if len(tokens) < 4:
                break
            fixid = tokens[1]
            group_name = tokens[2]
            delete_this_command = True
            assert(tokens[3].find('qmmm') == 0)
            if group_name in groups_needed:
                delete_this_command = False

            if delete_this_command:
                sys.stderr.write('WARNING: Ignoring line \n\"' +
                                 self.l_in_fix_qmmm[i_line].rstrip() + '\"\n')
                del self.l_in_fix_qmmm[i_line]
            else:
                self.l_in_fix_qmmm[i_line] = (' ' * self.indent) + ' '.join(tokens)
                i_line += 1



        sys.stderr.write(')\n\n')










    def Write(self, out_file):

        """
        ########################################
        ###  Now begin writing the template. ###
        ########################################
        """

        out_fname = None
        if isinstance(out_file, str):  # is "out_file" a file or a file name?
            out_fname = out_file
            out_file = open(out_fname, 'w')

        if not self.some_pair_coeffs_read:
            sys.stderr.write('Warning: No \"pair coeffs\" set.\n'
                             '         (No interactions between non-bonded atoms defined.)\n')
            self.no_warnings = False

        # sys.stderr.write('Writing ttree data to standard out.\n'
        #                 '       You can redirect this to a file using:\n'+
        #                 '   '+' '.join(sys.argv)+' > filename.ttree\n'
        #                 '        ----------------------\n')

        if self.mol_name != '':
            out_file.write(self.mol_name + ' {\n')

        if len(self.l_in_init) > 0:
            out_file.write('\n  ### LAMMPS commands for initialization\n'
                             '  ### (These can be overridden later.)\n\n')
            self.l_in_init.insert(0, (' ' * self.cindent) +
                             'write_once(\"' + in_init + '\") {')
            self.l_in_init.append((' ' * self.cindent) + '}\n')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_init))
        if len(self.l_in_settings) > 0:
            out_file.write('\n  ### LAMMPS commands for settings\n'
                             '  ### (These can be overridden later.)\n\n')
            self.l_in_settings.insert(0, (' ' * self.cindent) +
                                 'write_once(\"' + in_settings + '\") {')
            self.l_in_settings.append((' ' * self.cindent) + '}\n')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_settings))
            self.non_empty_output = True
        if len(self.l_in_masses) > 0:
            self.l_in_masses.insert(0, (' ' * self.cindent) +
                                    'write_once(\"' + in_settings + '\") {')
            self.l_in_masses.append((' ' * self.cindent) + '}\n')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_masses))
            self.non_empty_output = True


        if not self.ignore_coeffs:
            if self.remove_coeffs_from_data_file:
                if len(self.l_data_pair_coeffs) > 0:
                    for line in self.l_data_pair_coeffs:
                        #tokens = line.strip().split()
                        tokens = SplitQuotedString(line.strip(),
                                                   quotes='{',
                                                   delimiters=' \t\r\f\n',
                                                   escape='\\',
                                                   comment_char='#',
                                                   endquote='}')
                        atomtype_str = tokens[0]
                        self.l_in_pair_coeffs.append((' ' * self.cindent) +
                                                     '  pair_coeff ' +
                                                     atomtype_str + ' ' +
                                                     atomtype_str + ' ' +
                                                     ' '.join(tokens[1:]))
                    self.l_data_pair_coeffs = []
                if len(self.l_data_pairij_coeffs) > 0:
                    for line in self.l_data_pairij_coeffs:
                        #tokens = line.strip().split()
                        tokens = SplitQuotedString(line.strip(),
                                                   quotes='{',
                                                   delimiters=' \t\r\f\n',
                                                   escape='\\',
                                                   comment_char='#',
                                                   endquote='}')
                        atomtypeI_str = tokens[0]
                        atomtypeJ_str = tokens[1]
                        self.l_in_pair_coeffs.append((' ' * self.cindent) +
                                                     '  pair_coeff ' +
                                                     atomtypeI_str + ' ' +
                                                     atomtypeJ_str + ' ' +
                                                     ' '.join(tokens[2:]))
                        self.l_data_pairij_coeffs = []

            if len(self.l_in_pair_coeffs) > 0:
                self.l_in_pair_coeffs.insert(0, (' ' * self.cindent) +
                                             'write_once(\"' + in_settings + '\") {')
                self.l_in_pair_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_in_pair_coeffs) + '\n')
                self.non_empty_output = True

            if (self.remove_coeffs_from_data_file and (len(self.l_data_bond_coeffs) > 0)):
                for line in self.l_data_bond_coeffs:
                    self.l_in_bond_coeffs.append(
                        (' ' * self.cindent) + '  bond_coeff ' + line.strip())
                    self.l_data_bond_coeffs = []
            if len(self.l_in_bond_coeffs) > 0:
                self.l_in_bond_coeffs.insert(0, (' ' * self.cindent) +
                                             'write_once(\"' + in_settings + '\") {')
                self.l_in_bond_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_in_bond_coeffs) + '\n')
                self.non_empty_output = True

            if (self.remove_coeffs_from_data_file and (len(self.l_data_angle_coeffs) > 0)):
                for line in self.l_data_angle_coeffs:
                    self.l_in_angle_coeffs.append(
                        (' ' * self.cindent) + '  angle_coeff ' + line.strip())
                    self.l_data_angle_coeffs = []
                for line in self.l_data_bondbond_coeffs:
                    tokens = line.strip().split()
                    self.l_in_angle_coeffs.append(
                        (' ' * self.cindent) + '  angle_coeff ' + tokens[0] + ' bb ' + ' '.join(tokens[1:]))
                    self.l_data_bondbond_coeffs = []
                for line in self.l_data_bondangle_coeffs:
                    tokens = line.strip().split()
                    self.l_in_angle_coeffs.append(
                        (' ' * self.cindent) + '  angle_coeff ' + tokens[0] + ' ba ' + ' '.join(tokens[1:]))
                    self.l_data_bondangle_coeffs = []
            if len(self.l_in_angle_coeffs) > 0:
                self.l_in_angle_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + in_settings + '\") {')
                self.l_in_angle_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_in_angle_coeffs) + '\n')
                self.non_empty_output = True

            if (self.remove_coeffs_from_data_file and (len(self.l_data_dihedral_coeffs) > 0)):
                for line in self.l_data_dihedral_coeffs:
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + line.strip())
                    self.l_data_dihedral_coeffs = []

                for line in self.l_data_middlebondtorsion_coeffs:
                    tokens = line.strip().split()
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + tokens[0] + ' mbt ' + ' '.join(tokens[1:]))
                    self.l_data_middlebondtorsion_coeffs = []

                for line in self.l_data_endbondtorsion_coeffs:
                    tokens = line.strip().split()
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + tokens[0] + ' ebt ' + ' '.join(tokens[1:]))
                    self.l_data_endbondtorsion_coeffs = []

                for line in self.l_data_angletorsion_coeffs:
                    tokens = line.strip().split()
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + tokens[0] + ' at ' + ' '.join(tokens[1:]))
                    self.l_data_angletorsion_coeffs = []

                for line in self.l_data_angleangletorsion_coeffs:
                    tokens = line.strip().split()
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + tokens[0] + ' aat ' + ' '.join(tokens[1:]))
                    self.l_data_angleangletorsion_coeffs = []

                for line in self.l_data_bondbond13_coeffs:
                    tokens = line.strip().split()
                    self.l_in_dihedral_coeffs.append(
                        (' ' * self.cindent) + '  dihedral_coeff ' + tokens[0] + ' bb13 ' + ' '.join(tokens[1:]))
                    self.l_data_bondbond13_coeffs = []

            if len(self.l_in_dihedral_coeffs) > 0:
                self.l_in_dihedral_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + in_settings + '\") {')
                self.l_in_dihedral_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_in_dihedral_coeffs) + '\n')
                self.non_empty_output = True

            if (self.remove_coeffs_from_data_file and (len(self.l_data_improper_coeffs) > 0)):
                for line in self.l_data_improper_coeffs:
                    self.l_in_improper_coeffs.append(
                        (' ' * self.cindent) + '  improper_coeff ' + line.strip())
                    self.l_data_improper_coeffs = []

                for line in self.l_data_angleangle_coeffs:
                    tokens = line.strip().split()
                    self.l_in_improper_coeffs.append(
                        (' ' * self.cindent) + '  improper_coeff ' + tokens[0] + ' aa ' + ' '.join(tokens[1:]))
                    self.l_data_angleangle_coeffs = []

            if len(self.l_in_improper_coeffs) > 0:
                self.l_in_improper_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + in_settings + '\") {')
                self.l_in_improper_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_in_improper_coeffs) + '\n')
                self.non_empty_output = True


            # END OF "if not self.ignore_coeffs:"


        if self.non_empty_output:
            out_file.write('\n\n  ### DATA sections\n\n')

        if (len(self.l_data_masses) > 0) and (not self.ignore_masses):
            self.l_data_masses.insert(0, (' ' * self.cindent) +
                                 'write_once(\"' + data_masses + '\") {')
            self.l_data_masses.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_masses) + '\n')
            self.non_empty_output = True

        if not self.ignore_coeffs:
            if len(self.l_data_bond_coeffs) > 0:
                self.l_data_bond_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_bond_coeffs + '\") {')
                self.l_data_bond_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_bond_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_angle_coeffs) > 0:
                self.l_data_angle_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_angle_coeffs + '\") {')
                self.l_data_angle_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_angle_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_dihedral_coeffs) > 0:
                self.l_data_dihedral_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_dihedral_coeffs + '\") {')
                self.l_data_dihedral_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_dihedral_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_improper_coeffs) > 0:
                self.l_data_improper_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_improper_coeffs + '\") {')
                self.l_data_improper_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_improper_coeffs))
                out_file.write('\n')
                self.non_empty_output = True
            if len(self.l_data_pair_coeffs) > 0:
                self.l_data_pair_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_pair_coeffs + '\") {')
                self.l_data_pair_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_pair_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_pairij_coeffs) > 0:
                self.l_data_pairij_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_pairij_coeffs + '\") {')
                self.l_data_pairij_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_pairij_coeffs) + '\n')
                self.non_empty_output = True

            # class2 force fields:
            if len(self.l_data_bondbond_coeffs) > 0:
                self.l_data_bondbond_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_bondbond_coeffs + '\") {')
                self.l_data_bondbond_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_bondbond_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_bondangle_coeffs) > 0:
                self.l_data_bondangle_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_bondangle_coeffs + '\") {')
                self.l_data_bondangle_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_bondangle_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_middlebondtorsion_coeffs) > 0:
                self.l_data_middlebondtorsion_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_middlebondtorsion_coeffs + '\") {')
                self.l_data_middlebondtorsion_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_middlebondtorsion_coeffs) +'\n')
                self.non_empty_output = True
            if len(self.l_data_endbondtorsion_coeffs) > 0:
                self.l_data_endbondtorsion_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_endbondtorsion_coeffs + '\") {')
                self.l_data_endbondtorsion_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_endbondtorsion_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_angletorsion_coeffs) > 0:
                self.l_data_angletorsion_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_angletorsion_coeffs + '\") {')
                self.l_data_angletorsion_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_angletorsion_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_angleangletorsion_coeffs) > 0:
                self.l_data_angleangletorsion_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_angleangletorsion_coeffs + '\") {')
                self.l_data_angleangletorsion_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_angleangletorsion_coeffs)+'\n')
                self.non_empty_output = True
            if len(self.l_data_bondbond13_coeffs) > 0:
                self.l_data_bondbond13_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_bondbond13_coeffs + '\") {')
                self.l_data_bondbond13_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_bondbond13_coeffs) + '\n')
                self.non_empty_output = True
            if len(self.l_data_angleangle_coeffs) > 0:
                self.l_data_angleangle_coeffs.insert(
                    0, (' ' * self.cindent) + 'write_once(\"' + data_angleangle_coeffs + '\") {')
                self.l_data_angleangle_coeffs.append((' ' * self.cindent) + '}')
                out_file.write('\n')
                out_file.write('\n'.join(self.l_data_angleangle_coeffs) + '\n')
                self.non_empty_output = True


            # END OF "if not self.ignore_coeffs:"


        # automatic generation of bonded interactions by type:
        if len(self.l_data_angles_by_type) > 0:
            self.l_data_angles_by_type.insert(
                0, (' ' * self.cindent) + 'write_once(\"' + data_angles_by_type + '\") {')
            self.l_data_angles_by_type.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_angles_by_type) + '\n')
            self.non_empty_output = True
        if len(self.l_data_dihedrals_by_type) > 0:
            self.l_data_dihedrals_by_type.insert(
                0, (' ' * self.cindent) + 'write_once(\"' + data_dihedrals_by_type + '\") {')
            self.l_data_dihedrals_by_type.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_dihedrals_by_type) + '\n')
            self.non_empty_output = True
        if len(self.l_data_impropers_by_type) > 0:
            self.l_data_impropers_by_type.insert(
                0, (' ' * self.cindent) + 'write_once(\"' + data_impropers_by_type + '\") {')
            self.l_data_impropers_by_type.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_impropers_by_type) + '\n')
            self.non_empty_output = True

        if len(self.l_data_atoms) > 0:
            self.l_data_atoms.insert(0, (' ' * self.cindent) +
                                'write(\"' + data_atoms + '\") {')
            self.l_data_atoms.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_atoms) + '\n')
            self.non_empty_output = True
        else:
            sys.stderr.write('Warning: missing \"Atoms\" section.\n'
                             '         (Did you include a LAMMPS data file in your argument list?)\n')
            self.no_warnings = False

        # non-point-like particles
        if len(self.l_data_ellipsoids) > 0:
            self.l_data_ellipsoids.insert(
                0, (' ' * self.cindent) + 'write(\"' + data_ellipsoids + '\") {')
            self.l_data_ellipsoids.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_ellipsoids) + '\n')
        if len(self.l_data_lines) > 0:
            self.l_data_lines.insert(0, (' ' * self.cindent) +
                                'write(\"' + data_lines + '\") {')
            self.l_data_lines.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_lines) + '\n')
        if len(self.l_data_triangles) > 0:
            self.l_data_triangles.insert(0, (' ' * self.cindent) +
                                    'write(\"' + data_triangles + '\") {')
            self.l_data_triangles.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_triangles) + '\n')

        # DO NOT WRITE OUT VELOCITY DATA
        # (Why: because it makes it difficult to combine this molecular template
        #  with molecule templates from other sources which lack velocity data.
        #  LAMMPS (and topotools) will crash if the number of entries in the
        #  Velocities section of a data file does not match the number of atoms.)
        # COMMENTING OUT:
        # if len(self.l_data_velocities) > 0:
        #    self.l_data_velocities.insert(0, (' '*self.cindent)+'write(\"'+data_velocities+'\") {')
        #    self.l_data_velocities.append((' '*self.cindent)+'}')
        #    out_file.write('\n')
        #    out_file.write('\n'.join(self.l_data_velocities) + '\n')
        if len(self.l_data_bonds) > 0:
            if self.ignore_bond_types:
                self.l_data_bonds.insert(0, (' ' * self.cindent) +
                                         'write(\"' + data_bond_list + '\") {')
            else:
                self.l_data_bonds.insert(0, (' ' * self.cindent) +
                                         'write(\"' + data_bonds + '\") {')
            self.l_data_bonds.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_bonds) + '\n')
            self.non_empty_output = True
        if len(self.l_data_angles) > 0:
            self.l_data_angles.insert(0, (' ' * self.cindent) +
                                 'write(\"' + data_angles + '\") {')
            self.l_data_angles.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_angles) + '\n')
            self.non_empty_output = True
        if len(self.l_data_dihedrals) > 0:
            self.l_data_dihedrals.insert(0, (' ' * self.cindent) +
                                    'write(\"' + data_dihedrals + '\") {')
            self.l_data_dihedrals.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_dihedrals) + '\n')
            self.non_empty_output = True
        if len(self.l_data_impropers) > 0:
            self.l_data_impropers.insert(0, (' ' * self.cindent) +
                                    'write(\"' + data_impropers + '\") {')
            self.l_data_impropers.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_impropers) + '\n')
            self.non_empty_output = True
        if len(self.l_data_cmap) > 0:
            out_file.write('\n')
            self.l_data_cmap.insert(0, (' ' * self.cindent) +
                                    'write(\"' + data_cmap + '\") {')
            self.l_data_cmap.insert(0, '\n')
            self.l_data_cmap.insert(0, (' ' * self.cindent) +
                                'category $cmap(1, 1)\n')
            self.l_data_cmap.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_data_cmap) + '\n')
            self.non_empty_output = True

        if len(self.l_in_group) > 0:
            self.no_warnings = False
            self.l_in_group.insert(0, (' ' * self.cindent) +
                              'write(\"' + in_settings + '\") {')
            self.l_in_group.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_group) + '\n')
            # sys.stderr.write('######################################################\n'
            #                 'WARNING: One or more \"group\" commands appear to refer to relevant atoms.\n'
            #                 '         Please check to make sure that the group(s) generated by\n'
            #                 '         '+g_program_name+' contain the correct atoms.  (-Andrew 2017-10)\n'
            #                 '######################################################\n')
            assert(self.non_empty_output)

        if len(self.l_in_set) > 0:
            self.l_in_set.insert(0, ((' ' * self.cindent) +
                                'write(\"' + in_settings + '\") {'))
            self.l_in_set.append((' ' * self.cindent) + '} # end of list of \"set\" commands\n')
            out_file.write('\n')
            out_file.write((' ' * self.cindent) + '# list of \"set\" commands:\n')
            out_file.write('\n'.join(self.l_in_set) + '\n')

        if len(self.l_in_set_static) > 0:
            self.l_in_set_static.insert(0, ((' ' * self.cindent) +
                                       'write_once(\"' + in_settings + '\") {'))
            self.l_in_set_static.append((' ' * self.cindent) + '} # end of list of (static) \"set\" commands\n')
            out_file.write('\n')
            out_file.write((' ' * self.cindent) + '# list of (static) \"set\" commands:\n')
            out_file.write('\n'.join(self.l_in_set_static) + '\n')

        if len(self.l_in_fix_rigid) > 0:
            self.no_warnings = False
            self.l_in_fix_rigid.insert(0, (' ' * self.cindent) +
                                  'write(\"' + in_settings + '\") {')
            self.l_in_fix_rigid.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_fix_rigid) + '\n')
            sys.stderr.write('WARNING: \"fix rigid\" style command(s) applied to selected atoms.\n'
                             '         Please make sure that the fix group(s) are defined correctly.\n'
                             '######################################################\n')
            assert(self.non_empty_output)

        if len(self.l_in_fix_shake) > 0:
            self.no_warnings = False
            self.l_in_fix_shake.insert(0, (' ' * self.cindent) +
                                  'write(\"' + in_settings + '\") {')
            self.l_in_fix_shake.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_fix_shake) + '\n')
            sys.stderr.write('WARNING: \"fix shake\" style command(s) applied to selected atoms.\n'
                             '         Please check to make sure that the fix group(s) are defined correctly,\n'

                             '         and also check that the atom, bond, and angle types are correct.\n'
                             '######################################################\n')
            assert(self.non_empty_output)

        if len(self.l_in_fix_poems) > 0:
            self.no_warnings = False
            self.l_in_fix_poems.insert(0, (' ' * self.cindent) +
                                  'write(\"' + in_settings + '\") {')
            self.l_in_fix_poems.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_fix_poems) + '\n')
            sys.stderr.write('WARNING: \"fix poems\" style command(s) applied to selected atoms.\n'
                             '         Please make sure that the fix group(s) are defined correctly.\n'
                             '######################################################\n')
            assert(self.non_empty_output)

        if len(self.l_in_fix_qeq) > 0:
            self.no_warnings = False
            self.l_in_fix_qeq.insert(0, (' ' * self.cindent) +
                                'write(\"' + in_settings + '\") {')
            self.l_in_fix_qeq.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_fix_qeq) + '\n')
            sys.stderr.write('WARNING: \"fix qeq\" style command(s) applied to selected atoms.\n'
                             '         Please make sure that the fix group(s) are defined correctly.\n'
                             '######################################################\n')
            assert(self.non_empty_output)

        if len(self.l_in_fix_qmmm) > 0:
            self.no_warnings = False
            self.l_in_fix_qmmm.insert(0, (' ' * self.cindent) +
                                 'write(\"' + in_settings + '\") {')
            self.l_in_fix_qmmm.append((' ' * self.cindent) + '}')
            out_file.write('\n')
            out_file.write('\n'.join(self.l_in_fix_qmmm) + '\n')
            sys.stderr.write('WARNING: \"fix qmmm\" style command(s) applied to selected atoms.\n'
                             '         Please make sure that the fix group(s) are defined correctly.\n'
                             '######################################################\n')
            assert(self.non_empty_output)

        if self.mol_name != '':
            out_file.write('\n} # end of \"' + self.mol_name + '\" type definition\n')

        # if self.non_empty_output and self.no_warnings:
        if self.non_empty_output:
            sys.stderr.write('\n'
                             '# WARNING: Exotic (many-body) pair-styles and pair-styles with\n'
                             '#          unusual syntax (such hbond/dreiding) are not understood\n'
                             '#          by ' + g_program_name +
                             ' (...although they are supported by moltemplate).\n'
                             '#          Please look over the resulting LT file and check for errors.\n'
                             '#          Convert any remaining atom, bond, angle, dihedral, or improper id\n'
                             '#          or type numbers to the corresponding $ or @-style counter variables.\n'
                             '#          Feel free to report any bugs you find. (-Andrew Jewett 2017-10)\n')

        if out_fname != None:
            out_file.close()



def main():

    sys.stderr.write(g_program_name + ' v' +
                     g_version_str + ' ' + g_date_str + '\n')

    try:

        # Read the command line arguments to figure out which files
        # we will need to read (and also how to modify the output format).

        ltemp = Ltemplify(sys.argv[1:])

        # Note: The arguments in sys.argv contains the list of files that
        #       the user wants us to read and convert into moltemplate format.
        #       This list of files is now in the ltemp.input_files data member.

        # Using these settings, convert these files into a single file.
        # Send the content that file to sys.stdout.

        ltemp.Convert(sys.stdout,
                      ltemp.input_data_file,
                      ltemp.input_script_files)
                      

    except (InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return



if __name__ == '__main__':
    main()


