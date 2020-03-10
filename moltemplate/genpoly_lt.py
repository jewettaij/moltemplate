#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

"""
   Generate a moltemplate (.lt) file containing a definition of a polymer
   molecule whose monomers are located at the positions specified in
   "coords.raw" (a 3-column text file).  Monomers will be rotated so
   that they point in the direction connecting neighbors (r[i+1]-r[i])
   The user can specify the subunits to use when building the polymer,
   the atoms to to build bonds (and angles, and dihedrals) between monomers
   and the helical pitch of the polymer.  The output of this program is
   a text file in moltemplate (.lt) format containing the sequence of
   moltemplate commands needed to build this polymer molecule(s).  (One must
   then run moltemplate on this file to build the LAMMPS simulation files.)
       Multiple Polymers:
   To make it easier to create polymer melts, multiple polymers can be created
   from coordinates in the same file by using the "-cuts" command line argument.
      Encapsulation:
   If the "-polymer-name PolyName" command line option is given, then these
   moltemplate commands will be nested within the definition of a moltemplate
   object (named "PolyName", in this example. Later in your moltemplate files,
   you must remember to instantiate a copy of this moltemplate object using
   a command like "polymer = new PolyName"  Atoms within this object will
   share the same molecule-ID number.)  If multiple polymers are requested, then
   each of them will have their own polymer object.

"""


g_usage_msg = """
Usage:

   genpoly_lt.py  \\
      [-polymer-name pname] \\
      [-monomer-name mname] \\
      [-sequence sequence.txt] \\
      [-bond btype a1 a2] \\
      [-angle    atype a1 a2 a3 i1 i2 i3] \\
      [-dihedral dtype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-improper itype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-inherits ForceFieldObject] \\
      [-header "import \"monomer.lt\""] \\
      [-helix deltaphi] \\
      [-axis x,y,z] \\
      [-circular yes/no/connected] \\
      [-cuts cuts.txt] \\
      [-polymer-directions polarities.txt \\
      [-dir-indices ia ib] \\
      [-box paddingX,paddingY,paddingZ] \\
      < coords.raw > polymer.lt

"""


import sys
import random
from math import *


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



class GPSettings(object):

    def __init__(self):
        self.infile_name = ''
        self.direction_orig = [1.0, 0.0, 0.0]
        self.is_circular = False
        self.connect_ends = False
        self.delta_phi = 0.0
        #self.header = 'import \"forcefield.lt\"'
        self.header = ''
        self.name_monomer = 'Monomer'
        self.name_sequence_file = ''
        self.name_polymer = ''
        self.inherits = ''
        self.dir_index_offsets = (-1,1)
        self.cuts = []
        self.reverse_polymer_directions = []
        self.box_padding = None
        self.bonds_name = []
        self.bonds_type = []
        self.bonds_atoms = []
        self.bonds_index_offsets = []
        self.angles_name = []
        self.angles_type = []
        self.angles_atoms = []
        self.angles_index_offsets = []
        self.dihedrals_name = []
        self.dihedrals_type = []
        self.dihedrals_atoms = []
        self.dihedrals_index_offsets = []
        self.impropers_name = []
        self.impropers_type = []
        self.impropers_atoms = []
        self.impropers_index_offsets = []
        self.helix_angles_file = ''
        self.orientations_file = ''
        self.orientations_use_quats = False

    def ParseArgs(self, argv):
        i = 1
        while i < len(argv):
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

            if ((argv[i].lower() == '-in') or
                (argv[i].lower() == '-i')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name.\n')
                self.infile_name = argv[i + 1]
                del(argv[i:i + 2])
            elif argv[i].lower() == '-bond':
                if i + 3 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 4 strings.\n')
                # self.bonds_name.append(argv[i+1])
                self.bonds_type.append(argv[i + 1])
                self.bonds_atoms.append((argv[i + 2],
                                         argv[i + 3]))
                self.bonds_index_offsets.append((0, 1))
                del(argv[i:i + 4])
            elif argv[i].lower() == '-angle':
                if i + 7 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 5 strings and 3 integers.\n')
                # self.angles_name.append(argv[i+1])
                self.angles_type.append(argv[i + 1])
                self.angles_atoms.append((argv[i + 2],
                                          argv[i + 3],
                                          argv[i + 4]))
                self.angles_index_offsets.append((int(argv[i + 5]),
                                                  int(argv[i + 6]),
                                                  int(argv[i + 7])))
                if ((self.angles_index_offsets[-1][0] < 0) or
                    (self.angles_index_offsets[-1][1] < 0) or
                    (self.angles_index_offsets[-1][2] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3) must be >= 0\n')
                del(argv[i:i + 8])
            elif argv[i].lower() == '-dihedral':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 6 strings and 4 integers.\n')
                # self.dihedrals_name.append(argv[i+1])
                self.dihedrals_type.append(argv[i + 1])
                self.dihedrals_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                self.dihedrals_index_offsets.append((int(argv[i + 6]),
                                                     int(argv[i + 7]),
                                                     int(argv[i + 8]),
                                                     int(argv[i + 9])))
                if ((self.dihedrals_index_offsets[-1][0] < 0) or
                    (self.dihedrals_index_offsets[-1][1] < 0) or
                    (self.dihedrals_index_offsets[-1][2] < 0) or
                    (self.dihedrals_index_offsets[-1][3] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3 i4) must be >= 0\n')
                del(argv[i:i + 10])
            elif argv[i].lower() == '-improper':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 6 strings and 4 integers.\n')
                # self.impropers_name.append(argv[i+1])
                self.impropers_type.append(argv[i + 1])
                self.impropers_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                self.impropers_index_offsets.append((int(argv[i + 6]),
                                                     int(argv[i + 7]),
                                                     int(argv[i + 8]),
                                                     int(argv[i + 9])))
                if ((self.impropers_index_offsets[-1][0] < 0) or
                    (self.impropers_index_offsets[-1][1] < 0) or
                    (self.impropers_index_offsets[-1][2] < 0) or
                    (self.impropers_index_offsets[-1][3] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3 i4) must be >= 0\n')
                del(argv[i:i + 10])
            elif (argv[i].lower() == '-monomer-name'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.name_monomer = argv[i + 1]
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-sequence'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                self.name_sequence_file = argv[i + 1]
                del(argv[i:i + 2])

            elif (argv[i].lower() == '-cuts'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                try:
                    f = open(argv[i + 1], "r")
                except IOError:
                    raise InputError(
                        'Error: file ' + argv[i + 1] + ' could not be opened for reading\n')
                for line_orig in f:
                    line = line_orig.strip()
                    ic = line.find('#')
                    if ic != -1:
                        line = line[:ic]
                    else:
                        line = line.strip()
                    if len(line) > 0:
                        try:
                            self.cuts.append(int(line))
                        except ValueError:
                            raise InputError(
                                'Error: file ' + argv[i + 1] + ' should contain only nonnegative integers.\n')
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-polymer-directions'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                try:
                    f = open(argv[i + 1], "r")
                except IOError:
                    raise InputError(
                        'Error: file ' + argv[i + 1] + ' could not be opened for reading\n')
                for line_orig in f:
                    line = line_orig.strip()
                    ic = line.find('#')
                    if ic != -1:
                        line = line[:ic]
                    else:
                        line = line.strip()
                    if len(line) > 0:
                        try:
                            entry = line
                            if ((entry == '1') or
                                (entry == '+1') or
                                (entry == 'true') or
                                (entry == 'increasing')):
                                self.reverse_polymer_directions.append(False)
                            else:
                                self.reverse_polymer_directions.append(True)
                        except ValueError:
                            raise InputError(
                                'Error: file ' + argv[i + 1] + ' should contain only \"+1\" or \"-1\" on each line.\n')
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-polymer-name'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.name_polymer = argv[i + 1]
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-inherits'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.inherits = argv[i + 1]
                if self.inherits.find('inherits ') == 0:
                    self.inherits = ' ' + self.inherits
                else:
                    self.inherits = ' inherits ' + self.inherits
                if self.name_polymer == '':
                    self.name_polymer = 'Polymer'  # supply a default name
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-header'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string (usually in quotes)\n')
                if self.header != '':
                    self.header += '\n'
                self.header += argv[i + 1]
                del(argv[i:i + 2])
            elif argv[i].lower() == '-axis':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed ' +
                                     'by 3 numbers separated by commas (no spaces)\n')
                self.direction_orig = list(map(float, argv[i + 1].split(',')))
                del(argv[i:i + 2])
            elif argv[i].lower() == '-circular':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                if argv[i + 1].lower() == 'yes':
                    self.connect_ends = True
                    self.is_circular = True
                elif argv[i + 1].lower() == 'connected':
                    self.connect_ends = True
                    self.is_circular = False
                elif argv[i + 1].lower() == 'no':
                    self.connect_ends = False
                    self.is_circular = False
                else:
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                del(argv[i:i + 2])
            elif argv[i].lower() == '-helix':
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a number (angle in degrees)\n')
                self.delta_phi = float(argv[i + 1])
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-dir-indices'):
                if i + 2 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by two integers\n')
                self.dir_index_offsets = (int(argv[i + 1]), int(argv[i + 2]))
                if self.dir_index_offsets[0] == self.dir_index_offsets[1]:
                    raise InputError(
                        'Error: The two numbers following ' + argv[i] + ' must not be equal.\n')
                del(argv[i:i + 3])
            elif (argv[i].lower() == '-box'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed ' +
                                     'by 3 numbers separated by commas (no spaces)\n')
                self.box_padding = list(map(float, argv[i + 1].split(',')))
                if len(self.box_padding) == 1:
                    self.box_padding = self.box_padding * 3
                del(argv[i:i + 2])
            elif (argv[i].lower() in ('-orientations','-orientations')):
                self.orientations_file = argv[i+1]
                self.orientations_use_quats = False
                del(argv[i:i + 2])
            elif (argv[i].lower() in ('-quaternion','-quaternions')):
                self.orientations_file = argv[i+1]
                self.orientations_use_quats = True
                del(argv[i:i + 2])
            elif (argv[i].lower() in ('-helix-angle','-helix-angles')):
                self.helix_angles_file = argv[i+1]
                del(argv[i:i + 2])

            # elif ((argv[i][0] == '-') and (__name__ == '__main__')):
            #
            #    raise InputError('Error('+g_program_name+'):\n'+\
            #        'Unrecogized command line argument \"'+argv[i]+\
            #        '\"\n\n'+\
            #        __doc__)
            else:
                i += 1

        if ((len(self.reverse_polymer_directions) != 0) and 
            (len(self.reverse_polymer_directions) != len(self.cuts) + 1)):
            raise InputError('Error: The number of entries in the file you supplied to "-polymer-directions"\n'
                             '       does not equal the number of polymers (which is either 1, or 1 + the\n'
                             '       number of entries in the file you supplied to "-cuts", if applicable)\n')


        for b in range(0, len(self.bonds_type)):
            if len(self.bonds_type) > 1:
                self.bonds_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.bonds_name.append('genpoly')
        for b in range(0, len(self.angles_type)):
            if len(self.angles_type) > 1:
                self.angles_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.angles_name.append('genpoly')
        for b in range(0, len(self.dihedrals_type)):
            if len(self.dihedrals_type) > 1:
                self.dihedrals_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.dihedrals_name.append('genpoly')
        for b in range(0, len(self.impropers_type)):
            if len(self.impropers_type) > 1:
                self.impropers_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.impropers_name.append('genpoly')



class WrapPeriodic(object):
    """ Wrap() calculates the remainder of i % N.
        It turns out to be convenient to do this multiple times and later
        query whether i/N != 0 in any of them once (by checking bounds_err).

    """
    bounds_err = False

    @classmethod
    def Wrap(obj, i, N):
        if i // N != 0:
            obj.bounds_err = True
        return i % N

    def WrapF(obj, x, L):
        i = floor(x / L)
        if i != 0:
            obj.bounds_err = True
        return x - i * L


class GenPoly(object):
    """
        Read coordinates from a file, and generate a list of \"new\" commands
        in moltemplate format with the position of each monomer located
        at these positions, oriented appropriately, with bonds (and angles,
        dihedrals, etc...) connecting successive monomers together.
        By default, only a single polymer is created.
        However this class can create multiple polymers of different lengths.
        The list of coordinates for each polymer are saved separately within
        the "self.coords_multi" member.

    """

    def __init__(self):
        self.settings = GPSettings()
        self.coords_multi = []  # a list-of-list-of-lists of numbers Nxnx3
        self.name_sequence_multi = []
        self.direction_vects = []
        self.box_bounds_min = [0.0, 0.0, 0.0]
        self.box_bounds_max = [0.0, 0.0, 0.0]
        self.N = 0
        self.orientations_multi = []
        self.helix_angles_multi = []


    def ParseArgs(self, argv):
        """ 
        parse the argument list
        """
        # The command above will remove arguments from argv which are
        # understood by GPSettings.ParseArgs(argv).  
        # The remaining arguments will be handled below.
        self.settings.ParseArgs(argv)



    def ReadCoords(self, infile, ncolumns=3):
        """
        Read x y z coordinates from a multi-column ASCII text file.
        Store the coordinates in the self.coords_multi[][][] list.
        (Note: Information stored in self.settings.cuts and 
               self.settings.reverse_polymer_directions 
               is used to split the coordinate list into multiple lists
               and to determine the order that coordinates are read.)
        This function can also be used to read other (multi-column) ASCII text
        files by overriding the default value for "ncolumns".  (For example,
        to read files with 4 numbers on each line, set ncolumns=4.)
        """

        filename = ''
        if isinstance(infile, str):
            filename = infile
            try:
                infile = open(filename, 'r')
            except IOError:
                raise InputError(
                    'Error: file ' + filename +
                    ' could not be opened for reading\n')

        coords_multi = []  # (do not confuse this with "self.coords_multi")
        coords = []

        lines = infile.readlines()
        for i in range(0, len(lines)):
            tokens = lines[i].strip().split()
            if (len(tokens) == ncolumns):
                coords.append(list(map(float, tokens)))

        self.N = len(coords)
        if self.N < 2:
            err_msg = 'Error: Coordinate file must have at least 2 positions'
            raise InputError(err_msg+'.\n')

        # Did the caller ask us to split the polymer into multiple polymers?
        if len(self.settings.cuts) > 0:
            if (self.settings.cuts[-1] < self.N+1):
                self.settings.cuts.append(self.N + 1)
                self.settings.cuts.sort()
            i = 0
            for j in self.settings.cuts:
                coords_multi.append(coords[i:j])
                i = j
        else:
            coords_multi.append(coords)

        # Did the caller ask us to reverse the direction of any polymers?
        for i in range(0, len(self.settings.reverse_polymer_directions)):
            assert(i < len(self.coords_multi))
            if self.settings.reverse_polymer_directions[i]:
                self.coords_multi[i].reverse()

        if filename != '':
            # Then we opened a new file with that name.  We should close it now.
            infile.close()

        return coords_multi



    def ReadSequence(self, infile):
        """
        Read a sequence of monomer type names from a file.
        This function is similar to ReadCoords().
        """
        name_sequence = []

        filename = ''
        if isinstance(infile, str):
            filename = infile
            try:
                infile = open(filename, 'r')
            except IOError:
                raise InputError(
                    'Error: file ' + filename +
                    ' could not be opened for reading\n')

        for line_orig in infile:
            line = line_orig.strip()
            ic = line.find('#')
            if ic != -1:
                line = line[:ic]
            else:
                line = line.strip()
                if len(line) > 0:
                    name_sequence.append(line)

        N = len(name_sequence)

        # Did the caller ask us to split the polymer into multiple polymers?
        if len(self.settings.cuts) > 0:
            if (self.settings.cuts[-1] < N+1):
                self.settings.cuts.append(N + 1)
                self.settings.cuts.sort()
            i = 0
            for j in self.settings.cuts:
                self.name_sequence_multi.append(name_sequence[i:j])
                i = j
        else:
            self.name_sequence_multi.append(name_sequence)

        # Did the caller ask us to reverse the direction of any polymers?
        for i in range(0, len(self.settings.reverse_polymer_directions)):
            if self.settings.reverse_polymer_directions[i]:
                self.name_sequence_multi[i].reverse()

        if filename != '':
            # Then we opened a new file with that name.  We should close it now.
            infile.close()



    def ChooseDirections(self, coords):
        """
        Calculate the direction each monomer subunit should be pointing at:

        """

        N = len(coords)
        self.direction_vects = [[0.0, 0.0, 0.0] for i in range(0, N + 1)]

        if self.settings.is_circular:
            for i in range(0, N):
                # By default, the direction that monomer "i" is pointing is
                # determined by the position of the monomers before and after it
                # (at index i-1, and i+1).  More generally, we allow the user
                # to choose what these offsets are ("dir_index_offsets[")
                ia = WrapPeriodic.Wrap(i + self.settings.dir_index_offsets[0],
                                       N)
                ib = WrapPeriodic.Wrap(i + self.settings.dir_index_offsets[1],
                                       N)
                for d in range(0, 3):
                    self.direction_vects[i][d] = coords[ib][d] - coords[ia][d]
                        
        else:
            for i in range(1, N - 1):
                for d in range(0, 3):
                    self.direction_vects[i][d] = coords[
                        i + self.settings.dir_index_offsets[1]][d] - coords[
                            i + self.settings.dir_index_offsets[0]][d]

            for d in range(0, 3):
                self.direction_vects[0][d] = coords[1][d] - coords[0][d]
                self.direction_vects[N-1][d] = coords[N-1][d] - coords[N-2][d]

        # Optional: normalize the direction vectors

        for i in range(0, N):
            direction_len = 0.0
            for d in range(0, 3):
                direction_len += (self.direction_vects[i][d])**2
            direction_len = sqrt(direction_len)
            for d in range(0, 3):
                self.direction_vects[i][d] /= direction_len

        # Special case:  self.direction_vects[-1] is the direction that the original monomer
        # in "monomer.lt" was pointing.  (By default, 1,0,0 <--> the "x"
        # direction)

        self.direction_vects[-1] = self.settings.direction_orig



    def WriteLTFile(self, outfile):
        """ Write an moltemplate (.lt) file containing the definition of
        this polymer object.  (If multiple polymer objects were requested by
        the user (using the -cuts argument), then their definitions will
        appear nested within this object, and each of them will be
        instantiated once when the parent object is instantiated.)

        """

        # make sure len(genpoly.orientations_multi) == len(genpoly.coords_multi)
        if len(self.orientations_multi) == 0:
            self.orientations_multi = [[] for entry in self.coords_multi]
        # make sure len(self.helix_angles_multi) == len(self.coords_multi)
        if len(self.helix_angles_multi) == 0:
            self.helix_angles_multi = [[] for entry in self.coords_multi]

        outfile.write(self.settings.header + "\n\n\n")

        if len(self.coords_multi) == 1:
            self.WritePolymer(outfile,
                              self.settings.name_polymer +
                              self.settings.inherits,
                              self.coords_multi[0],
                              self.name_sequence_multi[0],
                              self.orientations_multi[0],
                              self.helix_angles_multi[0])
        else:
            if self.settings.name_polymer != '':
                outfile.write(self.settings.name_polymer + " {\n\n")
            outfile.write('# Definitions of individual polymers to follow\n\n')
            for i in range(0, len(self.coords_multi)):
                # insert ".../" in front of the inherits string
                ih_str = self.settings.inherits
                ih = ih_str.find('inherits ')
                if ih != -1:
                    ih += len('inherits ')
                    ih_str = ih_str[0:ih] + '.../' + ih_str[ih:]
                self.WritePolymer(outfile,
                                  self.settings.name_polymer + '_sub' + str(i + 1) +
                                  ih_str,
                                  self.coords_multi[i],
                                  self.name_sequence_multi[i],
                                  self.orientations_multi[i],
                                  self.helix_angles_multi[i])
            outfile.write('\n\n'
                          '# Now instantiate all the polymers (once each)\n\n')

            for i in range(0, len(self.coords_multi)):
                outfile.write('polymers[' + str(i) + '] = new ' +
                              self.settings.name_polymer + '_sub' + str(i + 1) + '\n')

            if self.settings.name_polymer != '':
                outfile.write('\n\n'
                              '}  # ' + self.settings.name_polymer + '\n\n')

        if self.settings.box_padding != None:
            for i in range(0, len(self.coords_multi)):
                # calculate the box big enough to collectively enclose
                # all of the coordinates (even multiple coordinate sets)
                self.CalcBoxBoundaries(self.coords_multi[i])
            self.WriteBoxBoundaries(outfile)



    def WritePolymer(self,
                     outfile,
                     name_polymer,
                     coords,
                     names_monomers,
                     orientations=[],
                     helix_angles=[]):
        """ Write a single polymer object to a file.
            This function is invoked by WriteLTFile()

        """
        N = len(coords)
        self.ChooseDirections(coords)

        if name_polymer != '':
            outfile.write(name_polymer + ' {\n'
                          '\n\n\n'
                          'create_var {$mol}\n'
                          '# The line above forces all monomer subunits to share the same molecule-ID\n'
                          '# (Note: Setting the molecule-ID number is optional and is usually ignored.)\n\n\n\n')

        outfile.write("""
# ------------ List of Monomers: ------------
#
# (Note: move(), rot(), and rotvv() commands control the position
#  of each monomer.  (See the moltemplate manual for an explanation
#  of what they do.)  Commands enclosed in push() are cumulative
#  and remain in effect until removed by pop().)



"""
                      )

        outfile.write("push(move(0,0,0))\n")

        phi = 0.0

        for i in range(0, N):
            #im1 = i-1
            # if im1 < 0 or self.settings.connect_ends:
            #    if im1 < 0:
            #        im1 += N
            outfile.write("pop()\n")
            if len(orientations) > 0:
                assert(len(orientations) == N)
                if self.settings.orientations_use_quats:
                    assert(len(orientations[i]) == 4)
                    outfile.write("push(quat(" +
                                  str(orientations[i][0]) + "," +
                                  str(orientations[i][1]) + "," +
                                  str(orientations[i][2]) + "," +
                                  str(orientations[i][3]) + "))\n")
                else:
                    assert(len(orientations[i]) == 9)
                    outfile.write("push(matrix(" +
                                  str(orientations[i][0]) + "," +
                                  str(orientations[i][1]) + "," +
                                  str(orientations[i][2]) + "," +
                                  str(orientations[i][3]) + "," +
                                  str(orientations[i][4]) + "," +
                                  str(orientations[i][5]) + "," +
                                  str(orientations[i][6]) + "," +
                                  str(orientations[i][7]) + "," +
                                  str(orientations[i][8]) + "))\n")
            else:
                # Otherwise, if no orientations were explicitly specified, then
                # infer the orientation from the direction of the displacement.
                outfile.write("push(rotvv(" +
                              str(self.direction_vects[i - 1][0]) + "," +
                              str(self.direction_vects[i - 1][1]) + "," +
                              str(self.direction_vects[i - 1][2]) + "," +
                              str(self.direction_vects[i][0]) + "," +
                              str(self.direction_vects[i][1]) + "," +
                              str(self.direction_vects[i][2]) + "))\n")
            outfile.write("push(move(" +
                          str(coords[i][0]) + "," +
                          str(coords[i][1]) + "," +
                          str(coords[i][2]) + "))\n")

            outfile.write("mon[" + str(i) + "] = new " +
                          names_monomers[i])

            # If requested, apply additional rotations about the polymer axis
            phi = 0.0
            if len(helix_angles) > 0:
                phi += helix_angles[i]
            elif self.settings.delta_phi != 0.0:
                phi = i * self.settings.delta_phi
            # Recall that self.direction_vects[-1] =
            # self.settings.direction_orig  (usually 1,0,0)
            outfile.write(".rot(" + str(phi) +
                          "," + str(self.settings.direction_orig[0]) +
                          "," + str(self.settings.direction_orig[1]) +
                          "," + str(self.settings.direction_orig[2]) +
                          ")\n")
            if len(orientations) > 0:
                outfile.write("pop()\n")


        outfile.write("pop()\n")
        if len(orientations) == 0:
            outfile.write("pop()\n"*N)

        assert(len(self.settings.bonds_name) ==
               len(self.settings.bonds_type) ==
               len(self.settings.bonds_atoms) ==
               len(self.settings.bonds_index_offsets))
        if len(self.settings.bonds_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Bonds\") {\n")
        WrapPeriodic.bounds_err = False
        for i in range(0, N):
            test = False
            for b in range(0, len(self.settings.bonds_type)):
                I = i + self.settings.bonds_index_offsets[b][0]
                J = i + self.settings.bonds_index_offsets[b][1]
                I = WrapPeriodic.Wrap(I, N)
                J = WrapPeriodic.Wrap(J, N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write(
                    "  $bond:" + self.settings.bonds_name[b] + str(i + 1))
                if len(self.settings.bonds_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @bond:" + self.settings.bonds_type[b] + " $atom:mon[" + str(I) + "]/" + self.settings.bonds_atoms[
                              b][0] + " $atom:mon[" + str(J) + "]/" + self.settings.bonds_atoms[b][1] + "\n")
        if len(self.settings.bonds_type) > 0:
            outfile.write("}  # write(\"Data Bonds\") {...\n\n\n")

        assert(len(self.settings.angles_name) ==
               len(self.settings.angles_type) ==
               len(self.settings.angles_atoms) ==
               len(self.settings.angles_index_offsets))
        if len(self.settings.angles_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Angles\") {\n")
        for i in range(0, N):
            for b in range(0, len(self.settings.angles_type)):
                I = i + self.settings.angles_index_offsets[b][0]
                J = i + self.settings.angles_index_offsets[b][1]
                K = i + self.settings.angles_index_offsets[b][2]
                I = WrapPeriodic.Wrap(I, N)
                J = WrapPeriodic.Wrap(J, N)
                K = WrapPeriodic.Wrap(K, N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write(
                    "  $angle:" + self.settings.angles_name[b] + str(i + 1))
                if len(self.settings.angles_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @angle:" + self.settings.angles_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.angles_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.angles_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.angles_atoms[b][2] +
                              "\n")
        if len(self.settings.angles_type) > 0:
            outfile.write("}  # write(\"Data Angles\") {...\n\n\n")

        assert(len(self.settings.dihedrals_name) ==
               len(self.settings.dihedrals_type) ==
               len(self.settings.dihedrals_atoms) ==
               len(self.settings.dihedrals_index_offsets))
        if len(self.settings.dihedrals_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Dihedrals\") {\n")
        for i in range(0, N):
            for b in range(0, len(self.settings.dihedrals_type)):
                I = i + self.settings.dihedrals_index_offsets[b][0]
                J = i + self.settings.dihedrals_index_offsets[b][1]
                K = i + self.settings.dihedrals_index_offsets[b][2]
                L = i + self.settings.dihedrals_index_offsets[b][3]
                I = WrapPeriodic.Wrap(I, N)
                J = WrapPeriodic.Wrap(J, N)
                K = WrapPeriodic.Wrap(K, N)
                L = WrapPeriodic.Wrap(L, N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write("  $dihedral:" +
                              self.settings.dihedrals_name[b] + str(i + 1))
                if len(self.settings.dihedrals_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @dihedral:" + self.settings.dihedrals_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.dihedrals_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.dihedrals_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.dihedrals_atoms[b][2] +
                              " $atom:mon[" + str(L) + "]/" + self.settings.dihedrals_atoms[b][3] +
                              "\n")
        if len(self.settings.dihedrals_type) > 0:
            outfile.write("}  # write(\"Data Dihedrals\") {...\n\n\n")

        assert(len(self.settings.impropers_name) ==
               len(self.settings.impropers_type) ==
               len(self.settings.impropers_atoms) ==
               len(self.settings.impropers_index_offsets))
        if len(self.settings.impropers_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Impropers\") {\n")
        for i in range(0, N):
            for b in range(0, len(self.settings.impropers_type)):
                I = i + self.settings.impropers_index_offsets[b][0]
                J = i + self.settings.impropers_index_offsets[b][1]
                K = i + self.settings.impropers_index_offsets[b][2]
                L = i + self.settings.impropers_index_offsets[b][3]
                I = WrapPeriodic.Wrap(I, N)
                J = WrapPeriodic.Wrap(J, N)
                K = WrapPeriodic.Wrap(K, N)
                L = WrapPeriodic.Wrap(L, N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write("  $improper:" +
                              self.settings.impropers_name[b] + str(i + 1))
                if len(self.settings.impropers_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @improper:" + self.settings.impropers_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.impropers_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.impropers_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.impropers_atoms[b][2] +
                              " $atom:mon[" + str(L) + "]/" + self.settings.impropers_atoms[b][3] +
                              "\n")
        if len(self.settings.impropers_type) > 0:
            outfile.write("}  # write(\"Data Impropers\") {...\n\n\n")

        if name_polymer != '':
            outfile.write("}  # " + name_polymer + "\n\n\n\n")

    def CalcBoxBoundaries(self, coords):
        N = len(coords)
        for i in range(0, N):
            for d in range(0, 3):
                if not self.box_bounds_min:
                    assert(not self.box_bounds_max)
                    self.box_bounds_min = [xd for xd in coords[i]]
                    self.box_bounds_max = [xd for xd in coords[i]]
                else:
                    if coords[i][d] > self.box_bounds_max[d]:
                        self.box_bounds_max[d] = coords[i][d]
                    if coords[i][d] < self.box_bounds_min[d]:
                        self.box_bounds_min[d] = coords[i][d]

    def WriteBoxBoundaries(self, outfile):
        for d in range(0, 3):
            self.box_bounds_min[d] -= self.settings.box_padding[d]
            self.box_bounds_max[d] += self.settings.box_padding[d]
        outfile.write("\n# ---------------- simulation box -----------------\n"
            
                      "# Now define a box big enough to hold a polymer with this (initial) shape\n"
                      "# (The user can override this later on.  This is the default box size.)"
                      "\n\n"
                      "write_once(\"Data Boundary\") {\n"
                      + str(self.box_bounds_min[0]) + "  " +
                      str(self.box_bounds_max[0]) + " xlo xhi\n"
                      + str(self.box_bounds_min[1]) + "  " +
                      str(self.box_bounds_max[1]) + " ylo yhi\n"
                      + str(self.box_bounds_min[2]) + "  " +
                      str(self.box_bounds_max[2]) + " zlo zhi\n"
                      "}\n\n\n")


def main():
    try:
        g_program_name = __file__.split('/')[-1]
        g_version_str = '0.1.0'
        g_date_str = '2019-12-13'
        sys.stderr.write(g_program_name + ' v' +
                         g_version_str + ' ' + g_date_str + '\n')
        argv = [arg for arg in sys.argv]
        genpoly = GenPoly()
        genpoly.ParseArgs(argv)
        # Any remain arguments?
        if len(argv) > 1:
            raise InputError('Error(' + g_program_name + '):\n' +
                             'Unrecogized command line argument \"' + argv[1] +
                             '\"\n\n' +
                             g_usage_msg)

        if genpoly.settings.infile_name != '':
            infile = open(genpoly.settings.infile_name, 'r')
        else:
            infile = sys.stdin
        outfile = sys.stdout

        # Read the coordinates
        genpoly.coords_multi = genpoly.ReadCoords(infile)

        # Did the user ask us to read a custom sequence of monomer type names?
        name_sequence_file = None
        if isinstance(genpoly.settings.name_sequence_file, str):
            if genpoly.settings.name_sequence_file != '':
                name_sequence_file = open(genpoly.settings.name_sequence_file,'r')
        if name_sequence_file:
            # Note: This will fill the contents of genpoly.name_sequence_multi
            genpoly.ReadSequence(name_sequence_file)
        else:
            # Otherwise just fill genpoly.name_sequence_multi with
            #  repeated copies of genpoly.settings.name_monomer
            #   (...using this ugly two-dimensional list-of-lists comprehension)
            genpoly.name_sequence_multi = [[genpoly.settings.name_monomer
                                            for j in
                                            range(0, len(genpoly.coords_multi[i]))]
                                           for i in range(0, len(genpoly.coords_multi))]

        # Did the user specify a file with a list of orientations?
        orientations_file = None
        if isinstance(genpoly.settings.orientations_file, str):
            if genpoly.settings.orientations_file != '':
                orientations_file = open(genpoly.settings.orientations_file, 'r')
        else:
            orientations_file = genpoly.settings.orientations_file
        if orientations_file:
            if genpoly.settings.orientations_use_quats:
                genpoly.orientations_multi = genpoly.ReadCoords(orientations_file, 4)
            else:
                genpoly.orientations_multi = genpoly.ReadCoords(orientations_file, 9)

        # Did the user specify a file with a list of helix angles?
        helix_angles_file = None
        if isinstance(genpoly.settings.helix_angles_file, str):
            if genpoly.settings.helix_angles_file != '':
                helix_angles_file = open(genpoly.settings.helix_angles_file, 'r')
        helix_angles_file = genpoly.settings.helix_angles_file
        if helix_angles_file:
            genpoly.helix_angles_multi = genpoly.ReadCoords(genpoly.settings.helix_angles_file, 1)
            # Because I borrowed the ReadCoords() function to read the file,
            # each "angle" is ended up being a list containing only one element
            # (the angle).  Convert this list into the number stored there.
            for i in range(0, len(genpoly.helix_angles_multi)):
                for j in range(0, len(genpoly.helix_angles_multi[i])):
                    assert(len(genpoly.helix_angles_multi[i][j]) == 1)
                    genpoly.helix_angles_multi[i][j] = genpoly.helix_angles_multi[i][j][0]

        # Now, check for polymer and sequence length inconsistencies:
        if (len(genpoly.coords_multi) != len(genpoly.name_sequence_multi)):
            raise InputError('Error(' +
                             g_program_name + '):\n' +
                             '      The coordinate file and sequence file have different lengths.\n')
        if ((len(genpoly.orientations_multi) > 0) and
            (len(genpoly.orientations_multi) != len(genpoly.coords_multi))):
            raise InputError('Error(' +
                             g_program_name + '):\n' +
                             '      The coordinate file and orientations/quats file have different lengths.\n')
        if ((len(genpoly.helix_angles_multi) > 0) and
            (len(genpoly.helix_angles_multi) != len(genpoly.coords_multi))):
            raise InputError('Error(' +
                             g_program_name + '):\n' +
                             '      The coordinate file and helix_angles file have different lengths.\n')
        for i in range(0, len(genpoly.coords_multi)):
            if len(genpoly.name_sequence_multi[i]) > 0:
                if (len(genpoly.coords_multi[i]) !=
                    len(genpoly.name_sequence_multi[i])):
                    raise InputError('Error(' +
                                     g_program_name + '):\n' +
                                     '      The coordinate file and sequence file have different lengths.\n')
            if ((len(genpoly.orientations_multi) > 0) and
                (len(genpoly.orientations_multi[i]) > 0)):
                if (len(genpoly.coords_multi[i]) !=
                    len(genpoly.orientations_multi[i])):
                    raise InputError('Error(' +
                                     g_program_name + '):\n' +
                                     '      The coordinate file and orientations/quats file have different lengths.\n')
            if ((len(genpoly.helix_angles_multi) > 0) and
                (len(genpoly.helix_angles_multi[i]) > 0)):
                if (len(genpoly.coords_multi[i]) !=
                    len(genpoly.helix_angles_multi[i])):
                    raise InputError('Error(' +
                                     g_program_name + '):\n' +
                                     '      The coordinate file and helix_angles file have different lengths.\n')



        # Convert all of this information to moltemplate (LT) format:
        genpoly.WriteLTFile(outfile)

        # Now close the input file
        if genpoly.settings.infile_name != '':  # <--if not reading from stdin
            infile.close()

    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()
