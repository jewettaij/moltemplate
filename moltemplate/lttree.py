#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.moltemplate.org
#         http://www.chem.ucsb.edu/~sheagroup
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2013, Regents of the University of California
# All rights reserved.

"""
lttree.py

lttree.py is an extension of the generic ttree.py program.
This version can understand and manipulate ttree-style templates which
are specialized for storing molecule-specific data for use in LAMMPS.

The main difference between lttree.py and ttree.py is:
Unlike ttree.py, lttree.py understands rigid-body movement commands like
"rot()" and "move()" which allows it to reorient and move each copy
of a molecule to a new location.  (ttree.py just ignores these commands.
Consequently LAMMPS input file (fragments) created with ttree.py have
invalid (overlapping) atomic coordinates and must be modified or aguemted
later (by loading atomic coordinates from a PDB file or an XYZ file).
lttree.py understands the "Data Atoms" section of a LAMMPS
data file (in addition to the various "atom_styles" which effect it).

Additional LAMMPS-specific features may be added in the future.

"""

g_program_name = __file__.split('/')[-1]  # ='lttree.py'
g_date_str = '2022-1-11'
g_version_str = '0.80.3'


import sys
from collections import defaultdict
import pkg_resources

try:
    from .ttree import BasicUISettings, BasicUIParseArgs, EraseTemplateFiles, \
        StackableCommand, PopCommand, PopRightCommand, PopLeftCommand, \
        PushCommand, PushLeftCommand, PushRightCommand, ScopeCommand, \
        WriteVarBindingsFile, StaticObj, InstanceObj, \
        BasicUI, ScopeBegin, ScopeEnd, WriteFileCommand, Render
    from .ttree_lex import InputError, TextBlock, DeleteLinesWithBadVars, \
        TemplateLexer, TableFromTemplate, VarRef, TextBlock, ErrorLeader, \
        SplitQuotedString
    from .lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid, \
        ColNames2Coords, ColNames2Vects, \
        data_atoms, data_prefix, data_masses, \
        data_velocities, data_ellipsoids, data_triangles, data_lines, \
        data_pair_coeffs, data_bond_coeffs, data_angle_coeffs, \
        data_dihedral_coeffs, data_improper_coeffs, data_bondbond_coeffs, \
        data_bondangle_coeffs, data_middlebondtorsion_coeffs, \
        data_endbondtorsion_coeffs, data_angletorsion_coeffs, \
        data_angleangletorsion_coeffs, data_bondbond13_coeffs, \
        data_angleangle_coeffs, data_bonds_by_type, data_angles_by_type, \
        data_dihedrals_by_type, data_impropers_by_type, \
        data_bonds, data_bond_list, data_angles, data_dihedrals, data_impropers, \
        data_boundary, data_pbc, data_prefix_no_space, in_init, in_settings, \
        in_prefix
    from .ttree_matrix_stack import AffineTransform, MultiAffineStack, \
        LinTransform, Matrix2Quaternion, MultQuat
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree import *
    from ttree_lex import *
    from lttree_styles import *
    from ttree_matrix_stack import *



try:
    unicode
except NameError:
    # Python 3
    basestring = unicode = str


class LttreeSettings(BasicUISettings):

    def __init__(self,
                 user_bindings_x=None,
                 user_bindings=None,
                 order_method='by_command'):

        BasicUISettings.__init__(self,
                                 user_bindings_x,
                                 user_bindings,
                                 order_method)

        # The following new member data indicate which columns store
        # LAMMPS-specific information.
        # The next 6 members store keep track of the different columns
        # of the "Data Atoms" section of a LAMMPS data file:
        self.column_names = []  # <--A list of column names (optional)
        self.ii_coords = []  # <--A list of triplets of column indexes storing coordinate data
        self.ii_vects = []  # <--A list of triplets of column indexes storing directional data
        #   (such as dipole or ellipsoid orientations)
        self.i_atomid = None  # <--An integer indicating which column has the atomid
        self.i_atomtype = None  # <--An integer indicating which column has the atomtype
        self.i_molid = None  # <--An integer indicating which column has the molid, if applicable
        self.print_full_atom_type_name_in_masses = False # <--how to print atom type names in the "Masses" section of a DATA file?




def LttreeParseArgs(argv, settings, main=False, show_warnings=True):
    # By default, include force_fields provided with the package
    argv.extend(["-import-path",
                 pkg_resources.resource_filename(__name__, 'force_fields/')])

    BasicUIParseArgs(argv, settings)

    # Loop over the remaining arguments not processed yet.
    # These arguments are specific to the lttree.py program
    # and are not understood by ttree.py:
    i = 1
    while i < len(argv):
        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
        if ((argv[i].lower() == '-atomstyle') or
                (argv[i].lower() == '-atom-style') or
                (argv[i].lower() == '-atom_style')):
            if i + 1 >= len(argv):
                raise InputError('Error(' + g_program_name + '): The ' + argv[i] + ' flag should be followed by a LAMMPS\n'
                                 '       atom_style name (or single quoted string containing a space-separated\n'
                                 '       list of column names such as: atom-ID atom-type q x y z molecule-ID.)\n')
            settings.column_names = AtomStyle2ColNames(argv[i + 1])
            sys.stderr.write('\n    \"' + data_atoms + '\" column format:\n')
            sys.stderr.write(
                '    ' + (' '.join(settings.column_names)) + '\n\n')
            settings.ii_coords = ColNames2Coords(settings.column_names)
            settings.ii_vects = ColNames2Vects(settings.column_names)
            settings.i_atomid, settings.i_atomtype, settings.i_molid = ColNames2AidAtypeMolid(
                settings.column_names)
            del(argv[i:i + 2])
        elif (argv[i].lower() == '-icoord'):
            if i + 1 >= len(argv):
                raise InputError('Error: ' + argv[i] + ' flag should be followed by list of integers\n'
                                 '       corresponding to column numbers for coordinates in\n'
                                 '       the \"' + data_atoms + '\" section of a LAMMPS data file.\n')
            ilist = argv[i + 1].split()
            if (len(ilist) % 3) != 0:
                raise InputError('Error: ' + argv[i] + ' flag should be followed by list of integers.\n'
                                 '       This is usually a list of 3 integers, but it can contain more.\n'
                                 '       The number of cooridnate columns must be divisible by 3,\n'
                                 '       (even if the simulation is in 2 dimensions)\n')
            settings.iaffinevects = []
            for i in range(0, len(ilist) / 3):
                cols = [int(ilist[3 * i]) + 1,
                        int(ilist[3 * i + 1]) + 1,
                        int(ilist[3 * i + 2]) + 1]
                settings.iaffinevects.append(cols)
            del(argv[i:i + 2])
        elif (argv[i].lower() == '-ivect'):
            if i + 1 >= len(argv):
                raise InputError('Error: ' + argv[i] + ' flag should be followed by list of integers\n'
                                 '       corresponding to column numbers for direction vectors in\n'
                                 '       the \"' + data_atoms + '\" section of a LAMMPS data file.\n')
            ilist = argv[i + 1].split()
            if (len(ilist) % 3) != 0:
                raise InputError('Error: ' + argv[i] + ' flag should be followed by list of integers.\n'
                                 '       This is usually a list of 3 integers, but it can contain more.\n'
                                 '       The number of cooridnate columns must be divisible by 3,\n'
                                 '       (even if the simulation is in 2 dimensions)\n')
            settings.ivects = []
            for i in range(0, len(ilist) / 3):
                cols = [int(ilist[3 * i]) + 1,
                        int(ilist[3 * i + 1]) + 1,
                        int(ilist[3 * i + 2]) + 1]
                settings.ivects.append(cols)
            del(argv[i:i + 2])
        elif ((argv[i].lower() == '-iatomid') or
              (argv[i].lower() == '-iid') or
              (argv[i].lower() == '-iatom-id')):
            if ((i + 1 >= len(argv)) or (not str.isdigit(argv[i + 1]))):
                raise InputError('Error: ' + argv[i] + ' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"' +
                                 data_atoms + '\" section of a\n'
                                 '       LAMMPS data file contains the atom id number (typically 1).\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            i_atomid = int(argv[i + 1]) - 1
            del(argv[i:i + 2])
        elif ((argv[i].lower() == '-iatomtype') or
              (argv[i].lower() == '-itype') or
              (argv[i].lower() == '-iatom-type')):
            if ((i + 1 >= len(argv)) or (not str.isdigit(argv[i + 1]))):
                raise InputError('Error: ' + argv[i] + ' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"' +
                                 data_atoms + '\" section of a\n'
                                 '       LAMMPS data file contains the atom type.\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            i_atomtype = int(argv[i + 1]) - 1
            del(argv[i:i + 2])
        elif ((argv[i].lower() == '-imolid') or
              (argv[i].lower() == '-imol') or
              (argv[i].lower() == '-imol-id') or
              (argv[i].lower() == '-imoleculeid') or
              (argv[i].lower() == '-imolecule-id')):
            if ((i + 1 >= len(argv)) or (not str.isdigit(argv[i + 1]))):
                raise InputError('Error: ' + argv[i] + ' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"' +
                                 data_atoms + '\" section of a\n'
                                 '       LAMMPS data file contains the molecule id number.\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            i_molid = int(argv[i + 1]) - 1
            del(argv[i:i + 2])

        elif (argv[i].lower() == '-full-comment-names'):
            settings.print_full_atom_type_name_in_masses = True
            del(argv[i:i + 1])

        elif (argv[i].lower() == '-short-comment-names'):
            settings.print_full_atom_type_name_in_masses = False
            del(argv[i:i + 1])

        elif (argv[i].find('-') == 0) and main:
            # elif (__name__ == "__main__"):
            raise InputError('Error(' + g_program_name + '):\n'
                             'Unrecogized command line argument \"' + argv[i] + '\"\n')
        else:
            i += 1


    if main:

        # Instantiate the lexer we will be using.
        #  (The lexer's __init__() function requires an openned file.
        #   Assuming __name__ == "__main__", then the name of that file should
        #   be the last remaining (unprocessed) argument in the argument list.
        #   Otherwise, then name of that file will be determined later by the
        #   python script which imports this module, so we let them handle it.)

        if len(argv) == 1:
            raise InputError('Error: This program requires at least one argument\n'
                             '       the name of a file containing ttree template commands\n')
        elif len(argv) == 2:
            try:
                # Parse text from the file named argv[1]
                settings.lex.infile = argv[1]  
                settings.lex.instream = open(argv[1], 'r')
            except IOError:
                sys.stderr.write('Error: unable to open file\n'
                                 '       \"' + argv[1] + '\"\n'
                                 '       for reading.\n')
                sys.exit(1)
            del(argv[1:2])

        else:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"' + arg + '\"' for arg in argv[1:]]
            raise InputError('Syntax Error(' + g_program_name + '):\n\n'
                             '       Problem with argument list.\n'
                             '       The remaining arguments are:\n\n'
                             '         ' + (' '.join(problem_args)) + '\n\n'
                             '       (The actual problem may be earlier in the argument list.\n'
                             '       If these arguments are source files, then keep in mind\n'
                             '       that this program can not parse multiple source files.)\n'
                             '       Check the syntax of the entire argument list.\n')

    if len(settings.ii_coords) == 0 and show_warnings:
        sys.stderr.write('########################################################\n'
                         '##            WARNING: atom_style unspecified         ##\n'
                         '## --> \"' + data_atoms + '\" column data has an unknown format ##\n'
                         '##              Assuming atom_style = \"full\"          ##\n'
                         #                 '########################################################\n'
                         #                 '## To specify the \"'+data_atoms+'\" column format you can:      ##\n'
                         #                 '##   1) Use the -atomstyle \"STYLE\"  argument         ##\n'
                         #                 '##      where \"STYLE\" is a string indicating a LAMMPS ##\n'
                         #                 '##      atom_style, including hybrid styles.(Standard ##\n'
                         #                 '##      atom styles defined in 2011 are supported.)   ##\n'
                         #                 '##   2) Use the -atomstyle \"COL_LIST\"    argument    ##\n'
                         #                 '##      where \"COL_LIST" is a quoted list of strings  ##\n'
                         #                 '##      indicating the name of each column.           ##\n'
                         #                 '##      Names \"x\",\"y\",\"z\" are interpreted as          ##\n'
                         #                 '##      atomic coordinates. \"mux\",\"muy\",\"muz\"         ##\n'
                         #                 '##      are interpreted as direction vectors.             ##\n'
                         #                 '##   3) Use the -icoord \"cx cy cz...\" argument        ##\n'
                         #                 '##      where \"cx cy cz\" is a list of integers        ##\n'
                         #                 '##      indicating the column numbers for the x,y,z   ##\n'
                         #                 '##      coordinates of each atom.                     ##\n'
                         #                 '##   4) Use the -ivect \"cmux cmuy cmuz...\" argument   ##\n'
                         #                 '##      where \"cmux cmuy cmuz...\" is a list of        ##\n'
                         #                 '##      integers indicating the column numbers for    ##\n'
                         #                 '##      the vector that determines the direction of a ##\n'
                         #                 '##      dipole or ellipsoid (ie. a rotateable vector).##\n'
                         #                 '##      (More than one triplet can be specified. The  ##\n'
                         #                 '##       number of entries must be divisible by 3.)   ##\n'
                         '########################################################\n')

        # The default atom_style is "full"
        settings.column_names = AtomStyle2ColNames('full')
        settings.ii_coords = ColNames2Coords(settings.column_names)
        settings.ii_vects = ColNames2Vects(settings.column_names)
        settings.i_atomid, settings.i_atomtype, settings.i_molid = ColNames2AidAtypeMolid(
            settings.column_names)

    return


def TransformAtomText(text, matrix, settings):
    """ Apply transformations to the coordinates and other vector degrees
    of freedom stored in the \"Data Atoms\" section of a LAMMPS data file.
    This is the \"text\" argument.
    The \"matrix\" stores the aggregate sum of combined transformations
    to be applied.

    """

    #sys.stderr.write('matrix_stack.M = \n'+ MatToStr(matrix) + '\n')

    lines = text.split('\n')

    for i in range(0, len(lines)):
        line_orig = lines[i]
        ic = line_orig.find('#')
        if ic != -1:
            line = line_orig[:ic]
            comment = ' ' + line_orig[ic:].rstrip('\n')
        else:
            line = line_orig.rstrip('\n')
            comment = ''

        # Split the line into words (columns) using whitespace delimeters
        columns = SplitQuotedString(line,
                                    quotes='{',
                                    endquote='}')

        if len(columns) > 0:
            if len(columns) == len(settings.column_names) + 3:
                raise InputError('Error: lttree.py does not yet support integer unit-cell counters \n'
                                 '   within the \"' + data_atoms + '\" section of a LAMMPS data file.\n'
                                 '   Instead please add the appropriate offsets (these offsets\n'
                                 '   should be multiples of the cell size) to the atom coordinates\n'
                                 '   in the data file, and eliminate the extra columns. Then try again.\n'
                                 '   (If you get this message often, email me and I\'ll fix this limitation.)')
            if len(columns) < len(settings.column_names):
                raise InputError('Error: The number of columns in your data file does not\n'
                                 '       match the LAMMPS atom_style you selected.\n'
                                 '       Use the -atomstyle <style> command line argument.\n')
            x0 = [0.0, 0.0, 0.0]
            x = [0.0, 0.0, 0.0]
            # Atomic coordinates transform using "affine" transformations
            # (translations plus rotations [or other linear transformations])
            for cxcycz in settings.ii_coords:
                for d in range(0, 3):
                    x0[d] = float(columns[cxcycz[d]])
                    AffineTransform(x, matrix, x0)  # x = matrix * x0 + b
                for d in range(0, 3):  # ("b" is part of "matrix")
                    columns[cxcycz[d]] = str(x[d])
            # Dipole moments and other direction-vectors
            # are not effected by translational movement
            for cxcycz in settings.ii_vects:
                for d in range(0, 3):
                    x0[d] = float(columns[cxcycz[d]])
                    LinTransform(x, matrix, x0)  # x = matrix * x0
                for d in range(0, 3):
                    columns[cxcycz[d]] = str(x[d])
        lines[i] = ' '.join(columns) + comment
    return '\n'.join(lines)



def TransformEllipsoidText(text, matrix, settings):
    """ Apply the transformation matrix to the quaternions represented
    by the last four numbers on each line.
    The \"matrix\" stores the aggregate sum of combined transformations
    to be applied and the rotational part of this matrix 
    must be converted to a quaternion.

    """

    #sys.stderr.write('matrix_stack.M = \n'+ MatToStr(matrix) + '\n')

    lines = text.split('\n')

    for i in range(0, len(lines)):
        line_orig = lines[i]
        ic = line_orig.find('#')
        if ic != -1:
            line = line_orig[:ic]
            comment = ' ' + line_orig[ic:].rstrip('\n')
        else:
            line = line_orig.rstrip('\n')
            comment = ''

        # Split the line into words (columns) using whitespace delimeters
        columns = SplitQuotedString(line,
                                    quotes='{',
                                    endquote='}')

        if len(columns) != 0:
            if len(columns) != 8:
                raise InputError('Error (lttree.py): Expected 7 numbers'
                                 + ' instead of '
                                 + str(len(columns))
                                 + '\nline:\n'
                                 + line
                                 + ' in each line of the ellipsoids\" section.\n"')
            q_orig = [float(columns[-4]),
                      float(columns[-3]),
                      float(columns[-2]),
                      float(columns[-1])]

            qRot = [0.0, 0.0, 0.0, 0.0]
            Matrix2Quaternion(matrix, qRot)

            q_new = [0.0, 0.0, 0.0, 0.0]
            MultQuat(q_new, qRot, q_orig)

            columns[-4] = str(q_new[0])
            columns[-3] = str(q_new[1])
            columns[-2] = str(q_new[2])
            columns[-1] = str(q_new[3])
            lines[i] = ' '.join(columns) + comment
    return '\n'.join(lines)



def CalcCM(text_Atoms,
           text_Masses=None,
           settings=None):
    types2masses = None
    # Loop through the "Masses" section: what is the mass of each atom type?
    if text_Masses != None:
        types2masses = {}
        lines = text_Masses.split('\n')
        for i in range(0, len(lines)):
            line = lines[i]
            # Split the line into words (columns) using whitespace delimeters
            columns = SplitQuotedString(line,
                                        quotes='{',
                                        endquote='}')
        if len(columns) == 2:
            atomtype = columns[0]
            m = float(columns[1])
            types2masses[atomtype] = m

    lines = text_Atoms.split('\n')
    # Pass 1 through the "Data Atoms" section: Determine each atom's mass
    if text_Masses != None:
        assert(settings != None)
        for i in range(0, len(lines)):
            line = lines[i]
            # Split the line into words (columns) using whitespace delimeters
            columns = SplitQuotedString(line,
                                        quotes='{',
                                        endquote='}')
            atomid = columns[settings.i_atomid]
            atomtype = columns[settings.i_atomtype]
            if atomtype not in types2masses[atomtype]:
                raise InputError('Error(lttree): You have neglected to define the mass of atom type: \"' + atomtype + '\"\n'
                                 'Did you specify the mass of every atom type using write(\"Masses\"){}?')
            atomid2mass[atomid] = atomtype2mass[atomtype]

    # Pass 2 through the "Data Atoms" section: Find the center of mass.
    for i in range(0, len(lines)):
        line = lines[i]
        # Split the line into words (columns) using whitespace delimeters
        columns = SplitQuotedString(line,
                                    quotes='{',
                                    endquote='}')
        if len(columns) > 0:
            if len(columns) == len(settings.column_names) + 3:
                raise InputError('Error: lttree.py does not yet support integer unit-cell counters (ix, iy, iz)\n'
                                 '   within the \"' + data_atoms + '\" section of a LAMMPS data file.\n'
                                 '   Instead please add the appropriate offsets (these offsets\n'
                                 '   should be multiples of the cell size) to the atom coordinates\n'
                                 '   in the data file, and eliminate the extra columns. Then try again.\n'
                                 '   (If you get this message often, email me and I\'ll fix this limitation.)')
            if len(columns) != len(settings.column_names):
                raise InputError('Error: The number of columns in your data file does not\n'
                                 '       match the LAMMPS atom_style you selected.\n'
                                 '       Use the -atomstyle <style> command line argument.\n')
            x = [0.0, 0.0, 0.0]
            if atomids2masses != None:
                m = atomids2masses[atomid]
            else:
                m = 1.0
            tot_m += m
            for cxcycz in settings.ii_coords:
                for d in range(0, 3):
                    x[d] = float(columns[cxcycz[d]])
                    tot_x[d] += x[d]
            # Note: dipole moments and other direction vectors don't effect
            #       the center of mass. So I commented out the loop below.
            # for cxcycz in settings.ii_vects:
            #    for d in range(0,3):
            #        v[d] = float(columns[cxcycz[d]])
        lines[i] = ' '.join(columns)

    xcm = [0.0, 0.0, 0.0]
    for d in range(0, 3):
        xcm[d] = tot_x[d] / tot_m
    return xcm





def AddAtomTypeComments(tmpl_list, substitute_vars, print_full_atom_type_names):
    """
    This ugly code attempts to parse the text in the "Masses" section
    of a LAMMPS DATA file, and append comments to the end of every line 
    defining the atom type.  Each comment will contain a string which stores 
    the name of the @atom-style variable (excluding the "@atom:" prefix).
    This is unfortunately complicated and messy because we have to do 
    this before we render the text.  (IE before we substutite numeric
    values into the variables.  Once we've rendered the text, 
    the variable names are discarded.)
    Therefore we have to work with a messy "tmpl_list" object
    which contains the text in a pre-rendered form.  The "tmpl_list" object
    is a list of alternating TextBlocks and VarRef objects.
    This function rebuilds this tmpl_list object, splitting it into separate
    lines (which it currently is not) and then adding comments to the end
    of each line (if there isn't one there already).  Finally it renders
    the resulting template and returns that text to the caller.
    """

    table = TableFromTemplate(tmpl_list,
                              [[' ', '\t', '\r'], '\n'],
                              [True, True])
    for i in range(0, len(table)):
        j = 0
        if isinstance(table[i][0], TextBlock):
            j += 1
        assert(hasattr(table[i], '__len__'))
        syntax_err = False
        if len(table[i]) == j+0:
            pass  # skip blank lines
        elif ((len(table[i]) > j+0) and
              isinstance(table[i][0], TextBlock) and
              (len(table[i][0].text) > 0) and
              (table[i][0].text == '#')):
            pass  # skip comment lines
        if ((len(table[i]) > j+1) and
            isinstance(table[i][j+0], VarRef) and
            isinstance(table[i][j+1], TextBlock)):
            var_ref = table[i][j+0]
            if print_full_atom_type_names:
                var_name = var_ref.prefix[0] + \
                            CanonicalDescrStr(var_ref.nptr.cat_name,
                                              var_ref.nptr.cat_node,
                                              var_ref.nptr.leaf_node,
                                              var_ref.srcloc)
            else:
                var_name = var_ref.nptr.leaf_node.name
            # remove the "@atom:" prefix before the variable name:
            if var_name.find('@atom:') == 0:
                var_name = var_name[6:]
            elif var_name.find('@/atom:') == 0:
                var_name = var_name[7:]
            new_comment = '  # ' + var_name
            if (len(table[i]) == j+2):
                table[i].append(TextBlock(new_comment,
                                          table[i][j+1].srcloc))
            else:
                assert(len(table[i]) > j+2)
                assert(isinstance(table[i][j+2], TextBlock))
                # If this line doesn't already contain a comment, then add one
                if table[i][j+2].text.find('#') == -1:
                    table[i][j+2].text += new_comment
                else:
                    # Insert a space between 2nd column and the comment
                    table[i][j+2].text = '  '+table[i][j+2].text

            # Also add spaces between any words within the comments.  This is
            # necessary because TableFromTemplate() removed all whitespace
            for k in range(j+3, len(table[i])):
                table[i][k].text = ' '+table[i][k].text
                
            # We must insert a space between the first and second columns
            # because TableFromTemplate() removes this whitespace separator.
            table[i].insert(j+1, TextBlock(' ', table[i][j+1].srcloc))
                    
        else:
            raise InputError('----------------------------------------------------\n' +
                             '     Syntax error near ' +
                             ErrorLeader(table[i][j+0].srcloc.infile,
                                         table[i][j+0].srcloc.lineno) + '\n'
                             '     The format is incorrect.\n')
        # Add a newline:
        table[i].append(TextBlock('\n',table[i][j+1].srcloc))

    # Now flatten the "table" (which is a list-of-lists) 
    # into a simple 1-dimensional list
    # (of alternating VarRefs and TextBlocks, in this case)

    templ_list = [entry for sublist in table for entry in sublist]

    # Note: This is equivalent to
    # templ_list = []
    # for sublist in table:
    #     for entry in sublist:
    #         templ_list.append(entry)
    # When building list comprehensions with multiple "for" tokens,
    # the outer loop comes first (ie "for sublist in table")

    # Now render this text and return it to the caller:
    return Render(templ_list, substitute_vars)



def _ExecCommands(command_list,
                  index,
                  global_files_content,
                  settings,
                  matrix_stack,
                  current_scope_id=None,
                  substitute_vars=True):
    """
    _ExecCommands():
    The argument "commands" is a nested list of lists of
    "Command" data structures (defined in ttree.py).

    Carry out the write() and write_once() commands (which
    write out the contents of the templates contain inside them).
    Instead of writing the files, save their contents in a string.

    The argument "global_files_content" should be of type defaultdict(list)
    It is an associative array whose key is a string (a filename)
    and whose value is a lists of strings (of rendered templates).

    """
    files_content = defaultdict(list)
    postprocessing_commands = []

    while index < len(command_list):
        command = command_list[index]
        index += 1

        # For debugging only
        if ((not isinstance(command, StackableCommand)) and
                (not isinstance(command, ScopeCommand)) and
                (not isinstance(command, WriteFileCommand))):
            sys.stderr.write(str(command) + '\n')

        if isinstance(command, PopCommand):
            assert(current_scope_id != None)
            if command.context_node == None:
                command.context_node = current_scope_id
            if isinstance(command, PopRightCommand):
                matrix_stack.PopRight(which_stack=command.context_node)
            elif isinstance(command, PopLeftCommand):
                matrix_stack.PopLeft(which_stack=command.context_node)
            else:
                assert(False)

        elif isinstance(command, PushCommand):
            assert(current_scope_id != None)
            if command.context_node == None:
                command.context_node = current_scope_id
            # Some commands are post-processing commands, and must be
            # carried out AFTER all the text has been rendered.  For example
            # the "movecm(0,0,0)" waits until all of the coordinates have
            # been rendered, calculates the center-of-mass, and then applies
            # a translation moving the center of mass to the origin (0,0,0).
            # We need to figure out which of these commands need to be
            # postponed, and which commands can be carried out now.
            # ("now"=pushing transformation matrices onto the matrix stack).
            # UNFORTUNATELY POSTPONING SOME COMMANDS MAKES THE CODE UGLY
            transform_list = command.contents.split('.')
            transform_blocks = []
            i_post_process = -1
            #  Example:  Suppose:
            #command.contents = '.rot(30,0,0,1).movecm(0,0,0).rot(45,1,0,0).scalecm(2.0).move(-2,1,0)'
            #  then
            #transform_list = ['rot(30,0,0,1)', 'movecm(0,0,0)', 'rot(45,1,0,0)', 'scalecm(2.0)', 'move(-2,1,0)']
            # Note: the first command 'rot(30,0,0,1)' is carried out now.
            # The remaining commands are carried out during post-processing,
            # (when processing the "ScopeEnd" command.
            #
            # We break up the commands into "blocks" separated by center-
            # of-mass transformations ('movecm', 'rotcm', or 'scalecm')
            #
            # transform_blocks = ['.rot(30,0,0,1)',
            #                     '.movecm(0,0,0).rot(45,1,0,0)',
            #                     '.scalecm(2.0).move(-2,1,0)']

            i = 0
            while i < len(transform_list):
                transform_block = ''
                while i < len(transform_list):
                    transform = transform_list[i]
                    i += 1
                    if transform != '':
                        transform_block += '.' + transform
                        transform = transform.split('(')[0]
                        if ((transform == 'movecm') or
                            (transform == 'rotcm') or
                            (transform == 'scalecm')):
                            break
                transform_blocks.append(transform_block)

            if len(postprocessing_commands) == 0:
                # The first block (before movecm, rotcm, or scalecm)
                # can be executed now by modifying the matrix stack.
                if isinstance(command, PushRightCommand):
                    matrix_stack.PushCommandsRight(transform_blocks[0].strip('.'),
                                                   command.srcloc,
                                                   which_stack=command.context_node)
                elif isinstance(command, PushLeftCommand):
                    matrix_stack.PushCommandsLeft(transform_blocks[0].strip('.'),
                                                  command.srcloc,
                                                  which_stack=command.context_node)
                # Everything else must be saved for later.
                postprocessing_blocks = transform_blocks[1:]
            else:
                # If we already encountered a "movecm" "rotcm" or "scalecm"
                # then all of the command blocks must be handled during
                # postprocessing.
                postprocessing_blocks = transform_blocks

            for transform_block in postprocessing_blocks:
                assert(isinstance(block, basestring))
                if isinstance(command, PushRightCommand):
                    postprocessing_commands.append(PushRightCommand(transform_block,
                                                                    command.srcloc,
                                                                    command.context_node))
                elif isinstance(command, PushLeftCommand):
                    postprocessing_commands.append(PushLeftCommand(transform_block,
                                                                   command.srcloc,
                                                                   command.context_node))

        elif isinstance(command, WriteFileCommand):

            # --- Throw away lines containin references to deleted variables:---

            # First: To edit the content of a template,
            #        you need to make a deep local copy of it
            tmpl_list = []
            for entry in command.tmpl_list:
                if isinstance(entry, TextBlock):
                    tmpl_list.append(TextBlock(entry.text,
                                               entry.srcloc))  # , entry.srcloc_end))
                else:
                    tmpl_list.append(entry)

            #     Now throw away lines with deleted variables

            DeleteLinesWithBadVars(tmpl_list)

            # --- Now render the text ---
            text = Render(tmpl_list,
                          substitute_vars)

            # ---- Coordinates of the atoms, must be rotated
            # and translated after rendering.
            # In addition, other vectors (dipoles, ellipsoid orientations)
            # must be processed.
            # This requires us to re-parse the contents of this text
            # (after it has been rendered), and apply these transformations
            # before passing them on to the caller.
            if command.filename == data_atoms:
                text = TransformAtomText(text, matrix_stack.M, settings)
            elif command.filename == data_ellipsoids:
                text = TransformEllipsoidText(text, matrix_stack.M, settings)
            if command.filename == data_masses:
                text = AddAtomTypeComments(tmpl_list,
                                           substitute_vars,
                                           settings.print_full_atom_type_name_in_masses)
            files_content[command.filename].append(text)

        elif isinstance(command, ScopeBegin):

            if isinstance(command.node, InstanceObj):
                if ((command.node.children != None) and
                        (len(command.node.children) > 0)):
                    matrix_stack.PushStack(command.node)

            # "command_list" is a long list of commands.
            # ScopeBegin and ScopeEnd are (usually) used to demarcate/enclose
            # the commands which are issued for a single class or
            # class instance.  _ExecCommands() carries out the commands for
            # a single class/instance.  If we reach a ScopeBegin(),
            # then recursively process the commands belonging to the child.
            index = _ExecCommands(command_list,
                                  index,
                                  files_content,
                                  settings,
                                  matrix_stack,
                                  command.node,
                                  substitute_vars)

        elif isinstance(command, ScopeEnd):
            if data_atoms in files_content:
                for ppcommand in postprocessing_commands:
                    if data_masses in files_content:
                        xcm = CalcCM(files_content[data_atoms],
                                     files_content[data_masses],
                                     settings)
                    else:
                        xcm = CalcCM(files_content[data_atoms])
                    if isinstance(ppcommand, PushRightCommand):
                        matrix_stack.PushCommandsRight(ppcommand.contents,
                                                       ppcommand.srcloc,
                                                       xcm,
                                                       which_stack=command.context_node)
                    elif isinstance(ppcommand, PushLeftCommand):
                        matrix_stack.PushCommandsLeft(ppcommand.contents,
                                                      ppcommand.srcloc,
                                                      xcm,
                                                      which_stack=command.context_node)
                    files_content[data_atoms] = \
                        TransformAtomText(files_content[data_atoms],
                                          matrix_stack.M, settings)
                    files_content[data_ellipsoids] = \
                        TransformEllipsoidText(files_content[data_ellipsoids],
                                               matrix_stack.M, settings)

                for ppcommand in postprocessing_commands:
                    matrix_stack.Pop(which_stack=command.context_node)
                    #(same as PopRight())

            if isinstance(command.node, InstanceObj):
                if ((command.node.children != None) and
                        (len(command.node.children) > 0)):
                    matrix_stack.PopStack()

            # "ScopeEnd"  means we're done with this class/instance.
            break

        else:
            assert(False)
            # no other command types allowed at this point

    # After processing the commands in this list,
    # merge the templates with the callers template list
    for filename, tmpl_list in files_content.items():
        global_files_content[filename] += \
            files_content[filename]

    return index


def ExecCommands(commands,
                 files_content,
                 settings,
                 substitute_vars=True):

    matrix_stack = MultiAffineStack()

    index = _ExecCommands(commands,
                          0,
                          files_content,
                          settings,
                          matrix_stack,
                          None,
                          substitute_vars)
    assert(index == len(commands))


def WriteFiles(files_content, suffix='', write_to_stdout=True):
    for filename, str_list in files_content.items():
        if filename != None:
            out_file = None
            if filename == '':
                if write_to_stdout:
                    out_file = sys.stdout
            else:
                out_file = open(filename + suffix, 'a')
            if out_file != None:
                out_file.write(''.join(str_list))
                if filename != '':
                    out_file.close()

    return


def main():
    """
    This is is a "main module" wrapper for invoking lttree.py
    as a stand alone program.  This program:

    1)reads a ttree file,
    2)constructs a tree of class definitions (g_objectdefs)
    3)constructs a tree of instantiated class objects (g_objects),
    4)automatically assigns values to the variables,
    5)and carries out the "write" commands to write the templates a file(s).

    """

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name + ' v' +
                     g_version_str + ' ' + g_date_str + ' ')
    sys.stderr.write('\n(python version ' + str(sys.version) + ')\n')
    if sys.version < '2.6':
        raise InputError(
            'Error: Alas, you must upgrade to a newer version of python.')

    try:

        #settings = BasicUISettings()
        #BasicUIParseArgs(sys.argv, settings)
        settings = LttreeSettings()
        LttreeParseArgs([arg for arg in sys.argv],  #(deep copy of sys.argv)
                        settings, main=True, show_warnings=True)

        # Data structures to store the class definitionss and instances
        g_objectdefs = StaticObj('', None)  # The root of the static tree
        # has name '' (equivalent to '/')
        g_objects = InstanceObj('', None)  # The root of the instance tree
        # has name '' (equivalent to '/')

        # A list of commands to carry out
        g_static_commands = []
        g_instance_commands = []

        BasicUI(settings,
                g_objectdefs,
                g_objects,
                g_static_commands,
                g_instance_commands)

        # Interpret the the commands.  (These are typically write() or
        # write_once() commands, rendering templates into text.
        # This step also handles coordinate transformations and delete commands.
        # Coordinate transformations can be applied to the rendered text
        # as a post-processing step.

        sys.stderr.write(' done\nbuilding templates...')

        files_content = defaultdict(list)

        ExecCommands(g_static_commands,
                     files_content,
                     settings,
                     False)
        ExecCommands(g_instance_commands,
                     files_content,
                     settings,
                     False)

        # Finally: write the rendered text to actual files.

        # Erase the files that will be written to:
        sys.stderr.write(' done\nwriting templates...')
        EraseTemplateFiles(g_static_commands)
        EraseTemplateFiles(g_instance_commands)

        # Write the files as templates
        # (with the original variable names present)
        WriteFiles(files_content, suffix=".template", write_to_stdout=False)

        # Write the files with the variables substituted by values
        sys.stderr.write(' done\nbuilding and rendering templates...')
        files_content = defaultdict(list)
        ExecCommands(g_static_commands, files_content, settings, True)
        ExecCommands(g_instance_commands, files_content, settings, True)
        sys.stderr.write(' done\nwriting rendered templates...\n')
        WriteFiles(files_content)
        sys.stderr.write(' done\n')

        # Now write the variable bindings/assignments table.
        sys.stderr.write('writing \"ttree_assignments.txt\" file...')
        # <-- erase previous version.
        open('ttree_assignments.txt', 'w').close()
        WriteVarBindingsFile(g_objectdefs)
        WriteVarBindingsFile(g_objects)
        sys.stderr.write(' done\n')

    except (ValueError, InputError) as err:
        if isinstance(err, ValueError):
            sys.stderr.write('Error converting string to numeric format.\n'
                             '      This sometimes means you have forgotten to specify the atom style\n'
                             '      (using the \"-atomstyle\" command).  Alternatively it could indicate\n'
                             '      that the moltemplate file contains non-numeric text in one of the\n'
                             '      .move(), .rot(), .scale(), .matrix(), or .quat() commands. If neither of\n'
                             '      these scenarios apply, please report this bug. (jewett.aij at gmail.com)\n')
            sys.exit(-1)
        else:
            sys.stderr.write('\n\n' + str(err) + '\n')
            sys.exit(-1)

    return


if __name__ == '__main__':
    main()
