#!/usr/bin/env python

g_program_name = __file__.split('/')[-1]
g_version_str  = '0.3.6'
g_date_str     = '2020-7-18'

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
      [-locations-random num_mods seed] \\
      [-width mod_width] \\
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




def FindNearestAvailableSite(i, # target index. look for a position closest to i
                             width,  # number of needed consecutive vacant sites
                             occupancy,     # an array of True,False values
                             is_periodic):  # consider "wrap around" indexing?
    """
    Look for an interval containing "width" vacant sites in the occupancy array
    (an array of True or False values) whose start is nearest to location i.
    """
    N = len(occupancy)
    # Check and see if site i is available.  If not, skip to the next site,
    # (either before or after this site).
    if is_periodic:
        j_stop = N // 2
    else:
        j_stop = max(-width+N-i, i)
    occupied = True
    for j in range(0, j_stop):
        occupied = True
        for s in (-1, 1):
            if i+s*j < 0:
                continue
            if i+s*j >= N:
                continue
            # Is location "i+s*j" available for an object of size "width"?
            # If so, all sites from [i+s*j, i+s*j+width) must not be occupied.
            occupied = False
            for d in range(0, width):
                if not is_periodic:
                    if i+s*j+d < 0:
                        occupied = True
                        continue
                    if i+s*j+d >= N:
                        occupied = True
                        continue
                if occupancy[(i+s*j+d) % N]:
                    occupied = True
                    break
            if not occupied:
                break
        if not occupied:
            break
    if occupied:
        return -1
    else:
        return i+s*j




def DistributePeriodic(widths,       # width of each object (>0)
                       occupancy,   # already occupied sites in the lattice
                       is_periodic=False, # is the lattice periodic?
                       offset=None):     # shift placment of first object?
    """
    Nm = len(widths)
    N = len(occupancy)
    Choose n integers from the interval [0,...,N-1] which are as
    evenly-spaced as possible, subject to the following constraints:
    1) Each of the n integers represents the position of an object which
       occupies "width" sites in a lattice with sites numbered [0, ..., N-1].
    2) Objects must avoid sites on the lattice which are already occupied.
    3) OPTIONAL: By default, the first integer is chosen arbitrarily.
       However you can specify the position of the first integer by providing
       an "offset" argument (an integer >=0). It will be the first of the n
       integers generated (or the closest available site to that location).
    The function returns a list of the n chosen integers.
    """

    Nm = len(widths)
    max_width = 0        
    if Nm > 0:
        max_width = max(widths)
    N = len(occupancy)
    locations = [-1 for im in range(0, Nm)]

    if is_periodic:
        Nreduced = N
    else:
        Nreduced = N - max_width
    # the next if statement is probably unnecessary (nobody invokes it this way)
    if offset == None:
        if is_periodic:
            offset = 0
        else:
            offset = Nreduced / (2*Nm)
    for im in range(0, Nm):
        i = offset + (N*im) // Nm  # next location?
        # If we didn't have to worry about occupancy, then we would de done now.
        # However if it is occupied, we have to find nearby unnoccupied sites:
        J = FindNearestAvailableSite(i, widths[im], occupancy, is_periodic)
        if J == -1:
            raise InputError('Error('+g_program_name+
                             '): Not enough available sites.\n')
        else:
            locations[im] = J
            for d in range(0, widths[im]):
                assert(occupancy[(J+d) % N] == False)
                occupancy[(J+d) % N] = True

    for im in range(0, Nm):          # error check: make sure that we remembered
        assert(locations[im] != -1) # to specify all the entries in locations[]

    return locations




def _DistributeRandom(widths,           # width of each object (>0)
                      occupancy,        # already occupied sites
                      is_periodic=False, # is the lattice periodic?
                      rand_seed=None):  # specify the random seed
    """
    Generate random non-overlapping integers in a 1-D lattice, taking care to
    avoid previously occupied lattice sites. Each integer has width "widths[im]",
    meaning that it occupies "widths[im]" sites on the lattice.  (The algorithm
    places each integer by inserting random amounts of space between them,
    taking into consideration their widths and the total lattice size.) Overlaps
    between integers with eachother and previously occupied sites are avoided.
    When there are previously occupied lattice sites, the algorithm used here
    is not smart enough to guarantee that it will place the objects in a truly
    random way, or even succeed in placing them at all.
    Running time: O(N), where N is the number of sites in the lattice.
    """

    Nm = len(widths)
    if Nm == 0:
        return []
    N = len(occupancy)
    locations = [-1 for im in range(0, Nm)]

    # "Nreduced" is the number of available sites in the "reduced" lattice.
    # Putting objects of width 1 (lattice site) in the reduced lattice
    # gives you the same number of choices that you would have by putting
    # objects of variable width in the original lattice.
    # So we will place width 1 objects in the reduced lattice, randomize their
    # position, and then figure out where they would be in the original lattice
    # by inserting widths[im]-1 new lattice sites following each object placment.
    # (Unfortunately, by placing objects in the reduced lattice, it's not
    #  obvious where they end up in the original lattice.  So its difficult
    #  to take into consideration which sites in the original lattice previously
    #  occupied and not available.  We will have to correct for overlaps with
    #  previously occupied sites later.  This is a limitation of this approach.)

    sum_widths = 0
    for im in range(0, Nm):
        sum_widths += widths[im]

    Nreduced = N - (sum_widths - Nm)  # size of the reduced lattice

    if Nreduced < Nm:
        raise InputError('Error('+g_program_name+'): Not enough space.\n')
    occupancy_reduced = [ -1 for im in range(0, Nreduced)]
    for im in range(0, Nm):
        occupancy_reduced[im] = im
    if rand_seed != None:
        random.seed(rand_seed)
    random.shuffle(occupancy_reduced)
    offset = 0
    if is_periodic:
        # (complicated boring detail)  By definition, each modification
        # occupies "self.widths[im]" monomers in the polymer.
        # In principle, the modification could occupy sites on the
        # polymer which cross the boundary between the last monomer
        # and the first monomer.  To allow this to happen, assume this does
        # not happen (as we have done so far), and then cyclically shift
        # the entries.  (The shift amount should be a random integer from
        # 0, max(widths)-1)
        offset = random.randint(0, max(widths)-1)

    # Index variables
    # im  =  which integer are we generating (ie. which object are we locating)
    # Ir =  which position in the reduced size lattice are we considering?
    # I  =  which position in the full size lattice are we considering?

    i = 0
    for ir in range(0, Nreduced):
        im = occupancy_reduced[ir]
        if im != -1:
            # Then "i" is the target site (in the original lattice) for
            # the im'th object we want to place.  Figure out whether site "i"
            # is available.  If not, find the nearest available site.
            J = FindNearestAvailableSite(i+offset,
                                         widths[im],
                                         occupancy,
                                         is_periodic)
            if J == -1:
                return None   #packing was unsuccessful during this attempt
            locations[im] = J
            for d in range(0, widths[im]):
                assert(occupancy[(J+d) % N] == False)
                occupancy[(J+d) % N] = True
            i += widths[im]
        else:
            i += 1

    for im in range(0, Nm):         # error check: make sure that we remembered
        assert(locations[im] != -1) # to specify all the entries in locations[]

    return locations



def DistributeRandom(widths,           # the width of each object (>0)
                     occupancy,        # already occupied sites
                     is_periodic=False, # is the lattice periodic?
                     rand_seed=None,   # specify the random seed
                     num_attempts=20): # number of randomly generated attempts
    """
    Generate random non-overlapping integers in a 1-D lattice, taking care to
    avoid previously occupied lattice sites. Each integer has width "widths[im]",
    meaning that it occupies "widths[im]" sites on the lattice.  (The algorithm
    places each integer by inserting random amounts of space between them,
    taking into consideration their widths and the total lattice size.) Overlaps
    between integers with eachother and previously occupied sites are avoided.
    When there are previously occupied lattice sites, the algorithm used here
    is not smart enough to guarantee that it will place the objects in a truly
    random way, or even succeed in placing them all.  So this function will
    attempt random placements a certain number of times before giving up.
    However, for sparsely occupied lattices, the number of attempts before
    success is O(1), and each attempt requires O(N) time (N=size of lattice)
    resulting in a running time of O(N), in that case.
    """
    if rand_seed == None:
        rand_seed = random.randrange(sys.maxsize)

    for a in range(0, num_attempts):
        occupancy_cpy = [i for i in occupancy] # a fresh copy of occupancy array
        L = _DistributeRandom(widths,
                              occupancy_cpy,
                              is_periodic,
                              rand_seed + a)
        if L != None:
            for i in range(0, len(occupancy)):  # if successful, then
                occupancy[i] = occupancy_cpy[i] # copy back into occupancy array
            break
    if L == None:
        raise InputError('Error('+g_program_name+
                         '): Not enough space.\n'+
                         '      (Quit after '+str(num_attempts)+
                         ' packing attempts.)\n')
    return L




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
        self.index_range = 0  # (this should equal to the maximum index offset)
        self.locations = []
        self.nmods = 0
        self.widths = [1]
        self.occupancy = []
        self.write_locations_file = ''
        self.write_occupancy_file = ''
        self.gen_locations_method = ''
        self.rand_seed = 0
        self.rand_num_attempts = 50
        self.periodic_offset = 0
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
            for line in lines:
                if ((len(line) == 0) or (line[0] == '#')):
                    continue
                I = int(line)
                if ((I < 0) or (self.N <= I)):
                    raise InputError('Error: Expected a number from 0 to '+str(self.N-1)+' on each line\n'
                                     '       of file "'+mod_locations_filename+'"\n'
                                     '       (NOTE: Indexing begins at 0, not 1.)\n')
                self.locations.append(I)
            self.nmods = len(self.locations)

        except (ValueError, InputError) as err:
            raise InputError('Error: Unable to read file "'+mod_locations_filename+'"\n'
                             '       -or- file has the wrong format (integers (>=0) on separate lines.)\n'
                             'Details:\n'+
                             str(err)+'\n')
            sys.exit(-1)



    def ParseArgs(self, argv):
        """ 
        Parse the argument list to define features which are specific to
        the specific polymer model I am using here.

        """
        mod_locations_filename = ''
        occupied_monomers = []
        i = 1
        while i < len(argv):

            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

            if (argv[i].lower() in ('-loc', '-locations')):
                if i + 1 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by a file name.\n')
                mod_locations_filename = argv[i+1]
                del(argv[i:i + 2])

            elif argv[i].lower() == '-length':
                if i + 1 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by a positive integer.\n')
                self.N = int(argv[i+1])
                del(argv[i:i+2])

            elif (argv[i].lower() in ('-polymer-name','-polymername')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by a file name.\n')
                self.polymer_name = argv[i + 1]
                del(argv[i:i + 2])

            elif argv[i].lower() in ('-locations-periodic',
                                     '-locationsperiodic',
                                     '-loc-periodic',
                                     '-locperiodic'):
                if i+2 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by 2 integers.\n')
                self.nmods = int(argv[i+1])
                self.gen_locations_method = 'periodic'
                self.periodic_offset = int(argv[i+2])
                del(argv[i:i + 3])

            elif argv[i].lower() in ('-locations-random', '-locationsrandom',
                                     '-loc-rand', '-locrand'):
                if i+2 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by 2 integers (n,seed)\n')
                self.gen_locations_method = 'random'
                self.nmods=int(argv[i+1])
                self.rand_seed = int(argv[i+2])
                del(argv[i:i + 3])

            elif argv[i].lower() in '-locations-random-attempts':
                if i+1 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by an integer\n')
                self.rand_num_attempts = int(argv[i+1])
                del(argv[i:i + 2])

            elif argv[i].lower() in ('-width', '-widths'):
                if i+1 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by an integer or a file name.\n')
                try:
                    try:
                        f = open(argv[i+1], 'r')
                        lines = f.readlines()
                        self.widths = []
                        for line in lines:
                            if line.strip() != '':
                                self.widths.append(int(line))
                        f.close()
                    except (IOError, OSError) as e:
                        self.widths=[int(argv[i+1])]
                except (ValueError) as e:
                    raise InputError('Error: '+argv[i]+' argument must be followed by a number (or a file of numbers).\n')
                for width in self.widths:
                    if width < 1:
                        raise InputError('Error: The number(s) supplied to the '+argv[i]+' argument must be >= 1.\n')
                del(argv[i:i + 2])

            elif argv[i].lower() in '-write-locations':
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by a file name.\n')
                self.write_locations_file = argv[i + 1]
                del(argv[i:i + 2])

            elif argv[i].lower() in '-write-occupancy':
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by a file name.\n')
                self.write_occupancy_file = argv[i + 1]
                del(argv[i:i + 2])

            elif ((argv[i].lower() == '-in') or
                (argv[i].lower() == '-i')):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by a file name.\n')
                self.infile_name = argv[i + 1]
                del(argv[i:i + 2])

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
                        'Error: The ' + argv[i] + ' flag should be followed by at least '+str(3+3*natoms)+
                        'arguments .\n')
                natoms = int(argv[i+1])
                if i+3+3*natoms >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' ' + argv[i+1] +
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
                    if offset > self.index_range:
                        self.index_range = offset
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
                        'Error: The ' + argv[i] + ' flag should be followed by at least'+str(5+2*natoms)+
                        'arguments .\n')
                natoms = int(argv[i+1])
                if i+5+2*natoms >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' ' + str(natoms) +
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
                    if offset > self.index_range:
                        self.index_range = offset
                self.fix_nbody_index_offsets.append(offsets)
                self.fix_nbody_params.append(argv[i+7+2*natoms])
                del(argv[i:i+7+2*natoms+1])

            elif argv[i].lower() == '-bond':
                if i + 5 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by 3 strings and 2 integers.\n')
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
                if 1 > self.index_range:
                    self.index_range = 1
                del(argv[i:i + 6])

            elif argv[i].lower() == '-angle':
                if i + 7 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by 4 strings and 3 integers.\n')
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
                    if offset > self.index_range:
                        self.index_range = offset
                self.angles_index_offsets.append(offsets)
                del(argv[i:i + 8])

            elif argv[i].lower() == '-dihedral':
                if i + 9 >= len(argv):
                    raise InputError('Error: The ' + argv[i] +
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
                    if offset > self.index_range:
                        self.index_range = offset
                self.dihedrals_index_offsets.append(offsets)
                del(argv[i:i + 10])

            elif argv[i].lower() == '-improper':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: The ' + argv[i] + ' flag should be followed by 5 strings and 4 integers.\n')
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
                    if offset > self.index_range:
                        self.index_range = offset
                self.impropers_index_offsets.append(offsets)
                del(argv[i:i + 10])

            elif argv[i].lower() == '-circular':
                if i + 1 >= len(argv):
                    raise InputError('Error: The ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                if argv[i + 1].lower() == 'yes':
                    self.connect_ends = True
                elif argv[i + 1].lower() == 'connected':
                    self.connect_ends = True
                elif argv[i + 1].lower() == 'no':
                    self.connect_ends = False
                else:
                    raise InputError('Error: The ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                del(argv[i:i + 2])

            elif argv[i].lower() in ('-read-occupancy', '-occupancy'):
                if i+1 >= len(argv):
                    raise InputError('Error: The '+argv[i]+' flag should be followed by a file name.\n')
                try:
                    try:
                        f = open(argv[i+1], 'r')
                        lines = f.readlines()
                        occupied_monomers = []
                        for line in lines:
                            if line.strip() != '':
                                occupied_monomers.append(int(line))
                        f.close()
                    except (IOError, OSError) as e:
                        raise InputError('Error: Unable to open file "'+argv[i+1]+'" for reading.\n')
                except (ValueError) as e:
                    raise InputError('Error: The "'+argv[i+1]+'" text file should contain\n'
                                     '       a list of numbers.  (One number per line.)\n')
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
        self.occupancy = [ False for i in range(0, self.N) ]

        if mod_locations_filename != '':
            self.LoadModLocations(mod_locations_filename)

        if self.nmods == 0:
            #raise InputError('Error: You have not specified the modifications you want to make.\n')
            sys.stderr.write('WARNING: You have not specified the modifications you want to make.\n')

        # The self.widths member should be a list containing self.nmods integers
        if len(self.widths) == 1:   # If the user supplied only one number, then
            self.widths = self.nmods * self.widths # use it for all the numbers.
        elif len(self.widths) != self.nmods:
            raise InputError('Error: The number of entries in the -widths file does ('+str(len(self.widths))+')\n'
                             '       does not match the number of modifications you requested (using the\n'
                             '       "-locations-random" or "-locations-periodic" arguments).\n')

        for i in range(0, self.nmods):
            # make sure each "width" is at least as large as self.index_range
            # (index_range equals the range of index offsets specified in the
            #  -angle, -dihedral, -improper, -set-atoms, -fix-nbody arguments.)
            if self.widths[i] < self.index_range:
                self.widths[i] = self.index_range

        if len(occupied_monomers) > 0:
            for i in occupied_monomers:
                if not ((0 <= i) and (i < self.N)):
                    raise InputError('Error: Encountered an invalide number in the occupancy file ('+str(i)+').\n'
                                     '       The integers in the occupancy file must lie between 0 and N-1,\n'
                                     '       where "N" is the number of monomers in the polymer ('+str(self.N)+').\n')
                assert(len(self.occupancy) == self.N)
                self.occupancy[i] = True

        if mod_locations_filename != '':
            pass
        if self.gen_locations_method == 'periodic':
            self.locations = DistributePeriodic(self.widths,
                                                self.occupancy,
                                                self.connect_ends,
                                                self.periodic_offset)

        elif self.gen_locations_method == 'random':
            self.locations = DistributeRandom(self.widths,
                                              self.occupancy,
                                              self.connect_ends,
                                              self.rand_seed,
                                              self.rand_num_attempts)
        self.nmods = len(self.locations)

        # Now that we know where all of the modifications will go, make sure
        # we update the "occupancy" array. (We might have done this already.)
        for i in range(0, self.nmods):
            for j in range(0, self.widths[i]):
                self.occupancy[self.locations[i]+j % self.N] = True




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

        if self.settings.nmods == 0:
            return

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
                for im in range(0, self.settings.nmods):
                    i = self.settings.locations[im]
                    ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                    if WrapPeriodic.bounds_err:
                        WrapPeriodic.bounds_err = False
                        if not self.settings.connect_ends:
                            continue
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
                for im in range(0, self.settings.nmods):
                    i = self.settings.locations[im]
                    ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                    if WrapPeriodic.bounds_err:
                        WrapPeriodic.bounds_err = False
                        if not self.settings.connect_ends:
                            continue
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
            for im in range(0, self.settings.nmods):
                i = self.settings.locations[im]
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                for b in range(0, len(self.settings.bonds_type)):
                    I = i + self.settings.bonds_index_offsets[b][0]
                    J = i + self.settings.bonds_index_offsets[b][1]
                    I = WrapPeriodic.Wrap(I, self.settings.N)
                    J = WrapPeriodic.Wrap(J, self.settings.N)
                    if WrapPeriodic.bounds_err:
                        WrapPeriodic.bounds_err = False
                        if not self.settings.connect_ends:
                            continue
                    if len(self.settings.bonds_type) > 1:
                        outfile.write('    $bond:gpm_bond'+str(b+1)+'_'+str(i+1))
                    else:
                        outfile.write('    $bond:gpm_bond_'+str(i+1))
                    outfile.write(' @bond:' + self.settings.bonds_type[b] +
                                  ' $atom:mon[' + str(I) + ']/' + self.settings.bonds_atoms[b][0] +
                                  ' $atom:mon[' + str(J) + ']/' + self.settings.bonds_atoms[b][1] +
                                  '\n')
            outfile.write('  }  # write("Data Bonds")\n')

        if len(self.settings.angles_type) > 0:
            outfile.write('\n')
            outfile.write('  write("Data Angles") {\n')
            WrapPeriodic.bounds_err = False
            for im in range(0, self.settings.nmods):
                i = self.settings.locations[im]
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
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
                    if len(self.settings.angles_type) > 1:
                        outfile.write('    $angle:gpm_angle'+str(b+1)+'_'+str(i+1))
                    else:
                        outfile.write('    $angle:gpm_angle_'+str(i+1))
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
            for im in range(0, self.settings.nmods):
                i = self.settings.locations[im]
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
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
                    if len(self.settings.dihedrals_type) > 1:
                        outfile.write('    $dihedral:gpm_dihedral'+str(b+1)+'_'+str(i+1))
                    else:
                        outfile.write('    $dihedral:gpm_dihedral_'+str(i+1))
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
            for im in range(0, self.settings.nmods):
                i = self.settings.locations[im]
                ip1 = WrapPeriodic.Wrap(i+1, self.settings.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
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
                    if len(self.settings.impropers_type) > 1:
                        outfile.write('    $improper:gpm_improper'+str(b+1)+'_'+str(i+1))
                    else:
                        outfile.write('    $improper:gpm_improper_'+str(i+1))
                    outfile.write(' @improper:' + self.settings.impropers_type[b] +
                                  ' $atom:mon[' + str(I) + ']/' + self.settings.impropers_atoms[b][0] +
                                  ' $atom:mon[' + str(J) + ']/' + self.settings.impropers_atoms[b][1] +
                                  ' $atom:mon[' + str(K) + ']/' + self.settings.impropers_atoms[b][2] +
                                  ' $atom:mon[' + str(L) + ']/' + self.settings.impropers_atoms[b][3] +
                                  '\n')
            outfile.write('  }  # write("Data Impropers")  \n')

        if self.settings.polymer_name != '':
            outfile.write('\n')
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

        if gen_poly_mod.settings.write_locations_file != '':
            f = open(gen_poly_mod.settings.write_locations_file, 'w')
            for im in range(0, gen_poly_mod.settings.nmods):
                f.write(str(gen_poly_mod.settings.locations[im])+'\n')
            f.close()

        if gen_poly_mod.settings.write_occupancy_file != '':
            f = open(gen_poly_mod.settings.write_occupancy_file, 'w')
            for i in range(0, gen_poly_mod.settings.N):
                if gen_poly_mod.settings.occupancy[i]:
                    f.write(str(i)+'\n')
            f.close()

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
