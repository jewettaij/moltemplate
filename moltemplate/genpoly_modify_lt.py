#!/usr/bin/env python

g_program_name = __file__.split('/')[-1]
g_version_str  = '0.1.1'
g_date_str     = '2020-4-15'

g_usage_msg = """

   Modify a polymer in specific locations specified by the user.
   Atom types can be modified, bonded interactions can be added or modified,
   and additional restraint fixes (eg. "fix restrain", "fix twist")
   can be added at these locations.

Usage:

   genpoly_modify_lt.py \\
      [-polymer-name pname] \\
      [-length num_monomers] \\
      [-locations filename] \\
      [-locations-periodic num_mods offset] \\
      [-locations-random num_mods] \\
      [-mod-width mod_width] \\
      [-bond btype a1 a2 i1 i2] \\
      [-angle    atype a1 a2 a3 i1 i2 i3] \\
      [-dihedral dtype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-improper itype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-set-atoms M filename attribute a1 ... aM i1 ... iM A1 ... Am] \\
      [-fix-nbody N filename fixname fixID group keyword a1 ... aN i1 ... iN params] \\
      [-circular yes/no/connected] \\
      >> polymer.lt

"""


# This arguments are not currently supported, and probably never will be:
#      [-cuts cuts.txt] \\
#      [-polymer-directions polarities.txt] \\
# ... so use this program only on single polymers.  (Don't use "-cut".)


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




class WrapPeriodic(object):
    """ 
    Wrap() returns the remainder of i % N.  I wrote this pointless-looking
    function as a convenient one-liner.  I will calculate this many times
    and later query whether i/N != 0 in any of them (by checking bounds_err).

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



class GPModSettings(object):

    def __init__(self):
        # Custom arguments only used for this polymer type:
        self.polymer_name = ''
        self.N = 0    # (number of monomers in the polymer (including end caps))
        self.connect_ends = False # are first and last monomers bonded together?
        self.end_padding = 0  # (this should equal to the maximum index offset)
        self.nmods = 0
        self.mods_evenly_spaced = False
        self.periodic_offset = 0
        self.mod_width = 1
        self.is_mod_here = [] # where should we make these modifications?
        self.escale = 1.0
        self.escale_name = ''
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
        self.fix_nbody_natoms = []
        self.fix_nbody_filename = []
        self.fix_nbody_fixID = []
        self.fix_nbody_group = []
        self.fix_nbody_fixname = []
        self.fix_nbody_keyword = []
        self.fix_nbody_atoms = []
        self.fix_nbody_index_offsets = []
        self.fix_nbody_params = []
        self.setatoms_natoms = []
        self.setatoms_filename = []
        self.setatoms_attribute_name = []
        self.setatoms_atoms = []
        self.setatoms_index_offsets = []
        self.setatoms_attributes = []



    def LoadModLocations(self, mod_locations_filename):
        try:
            loc_file = open(mod_locations_filename, 'r')
            lines = loc_file.readlines()
            for i in range(0, len(lines)):
                line = lines[i].strip()
                if ((len(line) == 0) or (line[0] == '#')):
                    continue
                n = int(lines[i])
                if ((i < 0) or (self.N <= i)):
                    raise InputError('Error: Expected a number from 0 to '+str(N-1)+' on each line\n'
                                     '       of file "'+mod_locations_file+'"\n'
                                     '       (NOTE: Indexing begins at 0, not 1.)\n')
                self.is_mod_here[n] = True
        except (ValueError, InputError) as err:
            sys.stderr.write('Error: Unable to read file "'+mod_locations_filename+'"\n'
                             '       -or- file has the wrong format (integers (>=0) on separate lines.)\n'
                             'Details:\n'+
                             str(err)+'\n')
            sys.exit(-1)



    def ChooseRandomModLocations(self):
        # "Navailable" is the number of available "sites" on the polymer
        # where the modification could be located.
        if self.connect_ends:
            Navailable = (self.N - (self.nmods*(self.mod_width-1)))
        else:
            Navailable = (self.N - ((self.nmods-1)*(self.mod_width-1) +
                                    max(self.mod_width-1, self.end_padding)))
        if Navailable < self.nmods:
            raise InputError('Error: Too many mods added (or mods are too wide) for a polymer of length '+str(self.N)+'\n'
                             '       (Try reducing the -nmods or -mod-width parameters.)\n')
        is_site_occupied = [ False for i in range(0, Navailable)]
        for i in range(0, self.nmods):
            is_site_occupied[i] = True
        random.shuffle(is_site_occupied)

        # Now fill the self.is_mod_here[] array, adding padding when necessary
        self.is_mod_here = [ False for i in range(0, self.N) ]
        j = 0
        for i in range(0, Navailable):
            if is_site_occupied[i]:
                assert(j + (self.mod_width-1) < self.N)
                self.is_mod_here[j] = True
                j += self.mod_width
            else:
                j += 1

        if self.connect_ends:
            # (complicated boring detail)  By definition, each modification
            # occupies "self.mod_width" monomers in the polymer.
            # In principle, the modification could occupy sites on the
            # polymer which cross the boundary between the last monomer
            # and the first monomer.  To allow this to happen, assume this does
            # not happen (as we have done so far), and then cyclically shift
            # the entries.  (The shift amount should be a random integer from
            # 0, self.mod_width-1)
            iso_cpy = [x for x in self.is_mod_here]
            shift = random.randint(0, self.mod_width-1)
            for i in range(0, self.N):
                self.is_mod_here[i] = iso_cpy[(i+shift) % self.N]

        return self.nmods


    def ChoosePeriodicModLocations(self):
        if self.connect_ends:
            Navailable = self.N
        else:
            Navailable = self.N - self.end_padding
        if self.periodic_offset == None:
            self.periodic_offset = Navailable / (2*self.nmods)
            if self.connect_ends:
                self.periodic_offset = 0
        for i in range(0, Navailable):
            ip1 = WrapPeriodic.Wrap(i+1, self.N)
            if ( ((-1+i-self.periodic_offset)*self.nmods) // Navailable <
                 ((i-self.periodic_offset)*self.nmods) // Navailable ):
                self.is_mod_here[i] = True
                sys.stderr.write('Modification made beginning at mon['+str(i)+']\n')


    def ParseArgs(self, argv):
        """ 
        Parse the argument list to define features which are specific to
        the specific polymer model I am using here.

        """
        mod_locations_filename = ''
        pmod = 0
        i = 1
        while i < len(argv):

            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

            if (argv[i].lower() in ('-loc', '-locations')):
                if i + 1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a file name.\n')
                mod_locations_filename = argv[i+1]
                self.mods_evenly_spaced = False
                del(argv[i:i + 2])

            #elif argv[i].lower() == '-pmod':
            #    if i+1 >= len(argv):
            #        raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
            #    try:
            #        pmod=int(argv[i+1])
            #    except ValueError:
            #        pmod=float(argv[i+1])
            #    del(argv[i:i + 2])

            elif argv[i].lower() in ('-locations-periodic',
                                     '-locationsperiodic',
                                     '-loc-periodic',
                                     '-locperiodic'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.nmods = int(argv[i+1])
                self.mods_evenly_spaced = True
                self.periodic_offset = int(argv[i+2])
                del(argv[i:i + 3])

            elif argv[i].lower() in ('-locations-random', '-locationsrandom',
                                     '-loc-rand', '-locrand'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.mods_evenly_spaced = False
                self.nmods=int(argv[i+1])
                del(argv[i:i + 2])

            elif argv[i].lower() in ('-mod-width', '-modwidth'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by an integer.\n')
                self.mod_width=int(argv[i+1])
                if self.mod_width < 1:
                    raise InputError('Error: '+argv[i]+' flag should be followed by an integer > 1.\n')
                del(argv[i:i + 2])

            elif argv[i].lower() == '-length':
                if i + 1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a positive integer.\n')
                self.N = int(argv[i+1])
                del(argv[i:i+2])

            elif (argv[i].lower() in ('-polymer-name','-polymername')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name.\n')
                self.polymer_name = argv[i + 1]
                del(argv[i:i + 2])

            elif ((argv[i].lower() == '-in') or
                (argv[i].lower() == '-i')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name.\n')
                self.infile_name = argv[i + 1]
                del(argv[i:i + 2])

            #elif ((argv[i].lower() == '-escale') or
            #      (argv[i].lower() == '-kT')):
            #    if i+2 >= len(argv):
            #        raise InputError('Error: '+argv[i]+' flag should be followed by a number and a string\n')
            #    self.escale = float(argv[i+1])
            #    self.escale_name = argv[i+2]
            #    del(argv[i:i+3])

            elif (argv[i].lower() in ('-set-atoms', '-set-atom', '-set')):

                """
                Parse arguments such as:
                -set-atoms M filename attribute a1 ... aM i1 ... iM A1 ... Am
                Eg:
                -atom-type 4 system.in.types type r c2 c2 r 0 0 1 1 Motor Motor Motor Motor
                """
                natoms = 1
                if i+3+3*natoms >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by at least '+str(3+3*natoms)+
                        'arguments .\n')
                natoms = int(argv[i+1])
                if i+3+3*natoms >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' ' + argv[i+1] +
                        ' flag should be followed by '+str(3+3*natoms-1)+
                        ' arguments .\n')
                self.setatoms_natoms.append(natoms)
                self.setatoms_filename.append(argv[i+2])
                self.setatoms_attribute_name.append(argv[i+3])
                atoms = []
                for j in range(0, natoms):
                    atoms.append(argv[i+4+j])
                self.setatoms_atoms.append(atoms)
                offsets = []
                for j in range(0, natoms):
                    offset = int(argv[i+4+natoms+j])
                    offsets.append(offset)
                    if offset < 0:
                        raise InputError(
                            'Error: ' + argv[i] + ' offset indices must be >= 0\n')
                    if offset > self.end_padding:
                        self.end_padding = offset
                self.setatoms_index_offsets.append(offsets)
                attributes = []
                for j in range(0, natoms):
                    attributes.append(argv[i+4+2*natoms+j])
                self.setatoms_attributes.append(attributes)
                del(argv[i:i+4+3*natoms])

            elif argv[i].lower() == '-fix-nbody':
                """
                Parse arguments such as:
                -fix-nbody 4 "fix_twist_62kTperturn.in" fxTw all twist "torque" r c2 c2 r 0 0 1 1 "36.9620502 0"
                """
                natoms = 1
                if i+5+2*natoms >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by at least'+str(5+2*natoms)+
                        'arguments .\n')
                natoms = int(argv[i+1])
                if i+5+2*natoms >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' ' + str(natoms) +
                        ' flag should be followed by '+str(5+2*natoms)+
                        'arguments .\n')
                self.fix_nbody_natoms.append(natoms)
                self.fix_nbody_filename.append(argv[i+2])
                self.fix_nbody_fixID.append(argv[i+3])
                self.fix_nbody_group.append(argv[i+4])
                self.fix_nbody_fixname.append(argv[i+5])
                self.fix_nbody_keyword.append(argv[i+6])
                atoms = []
                for j in range(0, natoms):
                    atoms.append(argv[i+7+j])
                self.fix_nbody_atoms.append(atoms)
                offsets = []
                for j in range(0, natoms):
                    offset = int(argv[i+7+natoms+j])
                    if offset < 0:
                        raise InputError(
                            'Error: ' + argv[i] + ' offset indices must be >= 0\n')
                    offsets.append(offset)
                    if offset > self.end_padding:
                        self.end_padding = offset
                self.fix_nbody_index_offsets.append(offsets)
                self.fix_nbody_params.append(argv[i+7+2*natoms])
                del(argv[i:i+7+2*natoms+1])

            elif argv[i].lower() == '-bond':
                if i + 5 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 3 strings and 2 integers.\n')
                # self.bonds_name.append(argv[i+1])
                self.bonds_type.append(argv[i + 1])
                self.bonds_atoms.append((argv[i + 2],
                                         argv[i + 3]))
                offsets = (int(argv[i + 4]),
                           int(argv[i + 5]))
                if ((offsets[0] < 0) or
                    (offsets[1] < 0)):
                    raise InputError('Error: ' + argv[i] +
                                     ' indices (i1 i2) must be >= 0\n')
                self.bonds_index_offsets.append(offsets)
                if 1 > self.end_padding:
                    self.end_padding = 1
                del(argv[i:i + 6])

            elif argv[i].lower() == '-angle':
                if i + 7 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 4 strings and 3 integers.\n')
                # self.angles_name.append(argv[i+1])
                self.angles_type.append(argv[i + 1])
                self.angles_atoms.append((argv[i + 2],
                                          argv[i + 3],
                                          argv[i + 4]))
                offsets = (int(argv[i + 5]),
                           int(argv[i + 6]),
                           int(argv[i + 7]))
                if ((offsets[0] < 0) or
                    (offsets[1] < 0) or
                    (offsets[2] < 0)):
                    raise InputError('Error: ' + argv[i] +
                                     ' indices (i1 i2 i3) must be >= 0\n')
                for offset in offsets:
                    if offset > self.end_padding:
                        self.end_padding = offset
                self.angles_index_offsets.append(offsets)
                del(argv[i:i + 8])

            elif argv[i].lower() == '-dihedral':
                if i + 9 >= len(argv):
                    raise InputError('Error: ' + argv[i] +
                                     ' flag should be followed by 5 strings and 4 integers.\n')
                # self.dihedrals_name.append(argv[i+1])
                self.dihedrals_type.append(argv[i + 1])
                self.dihedrals_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                offsets = (int(argv[i + 6]),
                           int(argv[i + 7]),
                           int(argv[i + 8]),
                           int(argv[i + 9]))
                if ((offsets[0] < 0) or
                    (offsets[1] < 0) or
                    (offsets[2] < 0) or
                    (offsets[3] < 0)):
                    raise InputError('Error: ' + argv[i] +
                                     ' indices (i1 i2 i3 i4) must be >= 0\n')
                for offset in offsets:
                    if offset > self.end_padding:
                        self.end_padding = offset
                self.dihedrals_index_offsets.append(offsets)
                del(argv[i:i + 10])

            elif argv[i].lower() == '-improper':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 5 strings and 4 integers.\n')
                # self.impropers_name.append(argv[i+1])
                self.impropers_type.append(argv[i + 1])
                self.impropers_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                offsets = (int(argv[i + 6]),
                           int(argv[i + 7]),
                           int(argv[i + 8]),
                           int(argv[i + 9]))
                if ((offsets[0] < 0) or
                    (offsets[1] < 0) or
                    (offsets[2] < 0) or
                    (offsets[3] < 0)):
                    raise InputError('Error: ' + argv[i] +
                                     ' indices (i1 i2 i3 i4) must be >= 0\n')
                for offset in offsets:
                    if offset > self.end_padding:
                        self.end_padding = offset
                self.impropers_index_offsets.append(offsets)
                del(argv[i:i + 10])

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

            elif ((argv[i].lower() == '-help') or (argv[i].lower() == '--help') 
                  or
                  (argv[i].lower() == '-?') or (argv[i].lower() == '--?')):
                sys.stderr.write('\n'+g_usage_msg+'\n')
                exit(0)
            
            elif ((argv[i][0] == '-') and (__name__ == '__main__')):
            
                raise InputError('Error('+g_program_name+'):\n'+\
                    'Unrecogized command line argument \''+argv[i]+\
                    '\'\n\n'+\
                    g_usage_msg)
            else:
                 i += 1

        if self.N == 0:
            raise InputError('Error: You must specify the length of the polymer\n'
                             '       using the "-length N" argument.\n')
        self.is_mod_here = [ False for i in range(0, self.N) ]

        if mod_locations_filename != '':
            self.LoadModLocations(mod_locations_filename)
        elif self.mods_evenly_spaced:
            self.ChoosePeriodicModLocations()
        else:
            self.ChooseRandomModLocations()






class GenPolyMod(object):

    def __init__(self):
        self.settings = GPModSettings()

    def ParseArgs(self, argv):
        """ 
        parse the argument list
        """
        # The command above will remove arguments from argv which are
        # understood by GFNettings.ParseArgs(argv).  
        # The remaining arguments will be handled below.
        self.settings.ParseArgs(argv)

    def WriteLTFile(self, outfile):

        if self.settings.polymer_name != '':
            outfile.write(self.settings.polymer_name + ' {\n')
            outfile.write('\n'
                          '  # (Note: This will augment the definition of "'+
                          self.settings.polymer_name+'".\n'
                          '  #        It will not overwrite or erase it.)\n'
                          '\n')

        if len(self.settings.setatoms_natoms) > 0:
            outfile.write('\n'
                          '  ########### Modifications using the "set" command ##########\n')

            for b in range(0, len(self.settings.setatoms_filename)):
                outfile.write('\n'
                              '  write("'+self.settings.setatoms_filename[b]+
                              '") {\n')
                WrapPeriodic.bounds_err = False
                for i in range(0, self.settings.N):
                    ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                    if WrapPeriodic.bounds_err:
                        WrapPeriodic.bounds_err = False
                        if not self.settings.connect_ends:
                            continue
                    if self.settings.is_mod_here[i]:
                        natoms = self.settings.setatoms_natoms[b]
                        for n in range(0, natoms):
                            I = i + self.settings.setatoms_index_offsets[b][n]
                            I = WrapPeriodic.Wrap(I, self.settings.N)
                            outfile.write('    set atom $atom:mon[' + str(I) + ']/' + self.settings.setatoms_atoms[b][n] +
                                          ' ' + self.settings.setatoms_attribute_name[b])
                            attribute=self.settings.setatoms_attributes[b][n]
                            if ((self.settings.setatoms_attribute_name[b] == 'type') and
                                (attribute.find('@atom:') != 0)):
                                attribute = '@atom:' + attribute
                            elif ((self.settings.setatoms_attribute_name[b] == 'mol') and
                                (attribute.find('$mol:') != 0)):
                                attribute = '$mol:' + attribute
                            outfile.write(' ' + attribute + '\n')
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                outfile.write('  }  # set atom '+self.settings.setatoms_attribute_name[b]+' ...\n'
                              '\n')

        # We can define the fixes that exert extra forces on some of the atoms
        # in the polymer, as well as where (which atoms) do they act on.

        if len(self.settings.fix_nbody_filename) > 0:
            outfile.write('\n'
                          '  ########### Modifications using the "fix" command ##########\n'
                          '  # Add nbody interactions mediated by fixes such as "fix restraint" and\n'
                          '  # "fix twist". (These fixes add forces between specific particles).\n')
            for b in range(0, len(self.settings.fix_nbody_filename)):
                outfile.write('\n'
                              '  write("'+self.settings.fix_nbody_filename[b] +
                              '") {\n'
                              '    fix '+self.settings.fix_nbody_fixID[b] +
                              ' ' + self.settings.fix_nbody_group[b] + 
                              ' ' + self.settings.fix_nbody_fixname[b])
                WrapPeriodic.bounds_err = False
                for i in range(0, self.settings.N):
                    ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                    if WrapPeriodic.bounds_err:
                        WrapPeriodic.bounds_err = False
                        if not self.settings.connect_ends:
                            continue
                    if self.settings.is_mod_here[i]:
                        outfile.write(' '+self.settings.fix_nbody_keyword[b])
                        natoms = self.settings.fix_nbody_natoms[b]
                        for n in range(0, natoms):
                            I = i + self.settings.fix_nbody_index_offsets[b][n]
                            I = WrapPeriodic.Wrap(I, self.settings.N)
                            outfile.write(' $atom:mon[' + str(I) + ']/' + self.settings.fix_nbody_atoms[b][n])
                        outfile.write(' '+self.settings.fix_nbody_params[b])
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                outfile.write('\n'
                              '  }  # write("fix ...\n'
                              '\n')


        if ((len(self.settings.bonds_type) > 0) or
            (len(self.settings.angles_type) > 0) or
            (len(self.settings.dihedrals_type) > 0) or
            (len(self.settings.impropers_type) > 0)):

            outfile.write('\n'
                          '  ########### Overriding dihedrals, angles, and impropers ##########\n'
                          '\n'
                          '  # Override the dihedrals (and angles, and impropers)\n'
                          '  # at certain locations.  This feature is sometimes used to add additional\n'
                          '  # bonded interactions to specific atoms in the polymer, or override\n'
                          '  # existing such interactions.  This feature is sometimes used\n'
                          '  # to turn off the dihedral (angle, improper) interactions\n'
                          '  # at locations in the polymer (that would otherwise interfere with the\n'
                          '  # behavior of other constraint forces (or twist motors) that we want to add.\n'
                          '  # One could (for example) replace a dihedral interaction\n'
                          '  # that constrains a dihedral angle, to one which exerts\n'
                          '  # no forces, allowing that dihedral angle to spin freely.\n'
                          '  # (One could also simply delete the dihedral angle, but that\n'
                          '  #  option is more difficult to implement and undo later.)\n')

        if len(self.settings.bonds_type) > 0:
            outfile.write('\n')
            outfile.write('  write("Data Bonds") {\n')
            WrapPeriodic.bounds_err = False
            for i in range(0, self.settings.N):
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                if self.settings.is_mod_here[i]:

                    for b in range(0, len(self.settings.bonds_type)):
                        I = i + self.settings.bonds_index_offsets[b][0]
                        J = i + self.settings.bonds_index_offsets[b][1]
                        I = WrapPeriodic.Wrap(I, self.settings.N)
                        J = WrapPeriodic.Wrap(J, self.settings.N)
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                        outfile.write('    $bond:gpm' + str(i + 1))
                        if len(self.settings.bonds_type) > 1:
                            outfile.write('_' + str(b + 1))
                        outfile.write(' @bond:' + self.settings.bonds_type[b] +
                                      ' $atom:mon[' + str(I) + ']/' + self.settings.bonds_atoms[b][0] +
                                      ' $atom:mon[' + str(J) + ']/' + self.settings.bonds_atoms[b][1] +
                                      '\n')
            outfile.write('  }  # write("Data Bonds")\n')

        if len(self.settings.angles_type) > 0:
            outfile.write('\n')
            outfile.write('  write("Data Angles") {\n')
            WrapPeriodic.bounds_err = False
            for i in range(0, self.settings.N):
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                if self.settings.is_mod_here[i]:

                    for b in range(0, len(self.settings.angles_type)):
                        I = i + self.settings.angles_index_offsets[b][0]
                        J = i + self.settings.angles_index_offsets[b][1]
                        K = i + self.settings.angles_index_offsets[b][2]
                        I = WrapPeriodic.Wrap(I, self.settings.N)
                        J = WrapPeriodic.Wrap(J, self.settings.N)
                        K = WrapPeriodic.Wrap(K, self.settings.N)
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                        outfile.write('    $angle:gpm' + str(i + 1))
                        if len(self.settings.angles_type) > 1:
                            outfile.write('_' + str(b + 1))
                        outfile.write(' @angle:' + self.settings.angles_type[b] +
                                      ' $atom:mon[' + str(I) + ']/' + self.settings.angles_atoms[b][0] +
                                      ' $atom:mon[' + str(J) + ']/' + self.settings.angles_atoms[b][1] +
                                      ' $atom:mon[' + str(K) + ']/' + self.settings.angles_atoms[b][2] +
                                      '\n')
            outfile.write('  }  # write("Data Angles")\n')

        if len(self.settings.dihedrals_type) > 0:
            outfile.write('\n')
            outfile.write('  write("Data Dihedrals") {\n')
            WrapPeriodic.bounds_err = False
            for i in range(0, self.settings.N):
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                if self.settings.is_mod_here[i]:

                    for b in range(0, len(self.settings.dihedrals_type)):
                        I = i + self.settings.dihedrals_index_offsets[b][0]
                        J = i + self.settings.dihedrals_index_offsets[b][1]
                        K = i + self.settings.dihedrals_index_offsets[b][2]
                        L = i + self.settings.dihedrals_index_offsets[b][3]
                        I = WrapPeriodic.Wrap(I, self.settings.N)
                        J = WrapPeriodic.Wrap(J, self.settings.N)
                        K = WrapPeriodic.Wrap(K, self.settings.N)
                        L = WrapPeriodic.Wrap(L, self.settings.N)
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                        outfile.write('    $dihedral:gpm' + str(i + 1))
                        if len(self.settings.dihedrals_type) > 1:
                            outfile.write('_' + str(b + 1))
                        outfile.write(' @dihedral:' + self.settings.dihedrals_type[b] +
                                      ' $atom:mon[' + str(I) + ']/' + self.settings.dihedrals_atoms[b][0] +
                                      ' $atom:mon[' + str(J) + ']/' + self.settings.dihedrals_atoms[b][1] +
                                      ' $atom:mon[' + str(K) + ']/' + self.settings.dihedrals_atoms[b][2] +
                                      ' $atom:mon[' + str(L) + ']/' + self.settings.dihedrals_atoms[b][3] +
                                      '\n')
            outfile.write('  }  # write("Data Dihedrals")\n')

        if len(self.settings.impropers_type) > 0:
            outfile.write('\n')
            outfile.write('  write("Data Impropers") {\n')
            WrapPeriodic.bounds_err = False
            for i in range(0, self.settings.N):
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                if self.settings.is_mod_here[i]:

                    for b in range(0, len(self.settings.impropers_type)):
                        I = i + self.settings.impropers_index_offsets[b][0]
                        J = i + self.settings.impropers_index_offsets[b][1]
                        K = i + self.settings.impropers_index_offsets[b][2]
                        L = i + self.settings.impropers_index_offsets[b][3]
                        I = WrapPeriodic.Wrap(I, self.settings.N)
                        J = WrapPeriodic.Wrap(J, self.settings.N)
                        K = WrapPeriodic.Wrap(K, self.settings.N)
                        L = WrapPeriodic.Wrap(L, self.settings.N)
                        if WrapPeriodic.bounds_err:
                            WrapPeriodic.bounds_err = False
                            if not self.settings.connect_ends:
                                continue
                        outfile.write('    $improper:gpm' + str(i + 1))
                        if len(self.settings.impropers_type) > 1:
                            outfile.write('_' + str(b + 1))
                        outfile.write(' @improper:' + self.settings.impropers_type[b] +
                                      ' $atom:mon[' + str(I) + ']/' + self.settings.impropers_atoms[b][0] +
                                      ' $atom:mon[' + str(J) + ']/' + self.settings.impropers_atoms[b][1] +
                                      ' $atom:mon[' + str(K) + ']/' + self.settings.impropers_atoms[b][2] +
                                      ' $atom:mon[' + str(L) + ']/' + self.settings.impropers_atoms[b][3] +
                                      '\n')
            outfile.write('  }  # write("Data Impropers")  \n')

        outfile.write('\n')

        if self.settings.polymer_name != '':
            outfile.write('}  # ' + self.settings.polymer_name + '\n')
            outfile.write('\n')





def main():
    try:
        sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')
        argv = [arg for arg in sys.argv]
        """
        if gen_poly_mod.settings.infile_name != '':
            infile = open(gen_poly_mod.settings.infile_name, 'r')
        else:
            infile = sys.stdin
        """
        outfile = sys.stdout
        gen_poly_mod = GenPolyMod()
        gen_poly_mod.ParseArgs(argv)
        # Any remaining arguments?
        if len(argv) > 1:
            raise InputError('Error('+g_program_name+'):\n'+
                             'Unrecogized command line argument "'+argv[1]+
                             '"\n\n'+
                             g_usage_msg)

        # Convert all of this information to moltemplate (LT) format:
        gen_poly_mod.WriteLTFile(outfile)

        """            
        # Now close the input file
        if gen_poly_mod.settings.infile_name != '': #<-if not reading from stdin
            infile.close()
        """

    #except (ValueError, InputError) as err:
    except (InputError) as err:
        sys.stderr.write('\n' + str(err))
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()
