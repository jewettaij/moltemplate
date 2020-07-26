#!/usr/bin/env python

g_program_name = __file__.split('/')[-1]
g_version_str  = '0.0.7'
g_date_str     = '2020-4-02'

usage_msg = """
Usage (example):

gen_nonbonded_tables.py \\
      [-types "C?,R C1motor,* Wall/Wall:DNAForceField/C*" ] \\
      [-label "DNA_U0=1kT"] \\
      [-table-file "table_DNA_U0=1kT.dat"] \\
      [-script-file "pair_DNA_U0=1kT.dat"] \\
      [-pair-style "table"] \\
      [-barrier-height 1] \\
      [-escale 0.5961621] \\
      [-lj epsilon sigma Lcoeff] \\
      [-ljpq epsilon sigma Lcoeff p q] \\
      [-ljrmax rmax_ratio] \\
      [-rmin rmin] \\
      [-charge q] \\
      [-Ldebye 1.5] \\
      [-Ldebye-cut 7.5] \\

"""


import sys
import random
from math import *
import calc_table_nonbonded
import re        



class GNBTSettings(object):

    def __init__(self):
        self.atom_type_pairs = []
        self.table_filename = ''
        self.script_filename = ''
        self.pair_style_name = ''
        self.object_name = ''
        self.label = ''
        self.charge1 = 0.0
        self.charge2 = None
        self.Ldebye = 1.0 # Debye length
        self.Ldebyecut = 5.0
        self.rshiftQ = 0.0
        self.eps = 0.0
        self.sig = 1.0
        self.Lcoeff = 0.0
        self.Ntable = 64
        self.rmin = 0.0
        self.rmax_LJ = 0.0
        self.rmax_LJ_ratio = 1.0 # location of minima (purely repulsive)
        self.rshiftLJ = 0.0
        self.U0 = float('inf')
        self.escale = 1.0
        self.escale_name = ''
        #kB = 0.001987207 # (kCal/mole) / degreeK
        #Tkelvin  = 300.0
        #self.escale = kB*Tkelvin

    def ParseArgs(self, argv):
        """ 
        Parse the argument list to define features which are specific to
        the specific polymer model I am using here.
        (General settings for polymers are handled by GPSettings.ParseArgs())

        """
        i = 1
        while i < len(argv):
            #sys.stderr.write('argv['+str(i)+'] = "'+argv[i]+'"\n')
            if (argv[i].lower() in ('-types','-atom-types')):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a string\n')
                type_str = argv[i+1]
                # the next line works, but it does not ignore \ escaped spaces
                #type_pair = type_str.split(' ')
                type_list = re.split(r'(?<!\\) ', type_str) #considers \ escape
                if (len(type_list) % 2 != 0):
                    raise InputError('Error: Expected an even number of space-delimited entries in a quoted string\n'
                                     '       following the "'+argv[i]+'" argument.\n'
                                     'Example: '+argv[i]+' "C1 C2"\n')
                for j in range(0, len(type_list)):
                    # If the string r'\ ' appears anywhere, replace it with ' '
                    type_list[j] = type_list[j].replace(r'\ ', ' ')
                type_pairs = []
                for j in range(0, len(type_list)//2):
                    type_pairs.append((type_list[2*j], type_list[2*j+1]))
                self.atom_type_pairs += type_pairs
                del(argv[i:i+2])
            elif (argv[i].lower() == '-table-file'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.table_filename  = argv[i+1]
                del(argv[i:i+2])
            elif (argv[i].lower() == '-script-file'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.script_filename  = argv[i+1]
                del(argv[i:i+2])
            elif (argv[i].lower() == '-pair-style'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.pair_style_name  = argv[i+1]
                del(argv[i:i+2])
            elif (argv[i].lower() == '-label'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.label  = argv[i+1]
                del(argv[i:i+2])
            elif (argv[i].lower() in ('-object-name','-object')):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.object_name  = argv[i+1]
                del(argv[i:i+2])
            elif (argv[i].lower() == '-rmin'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.rmin = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-Ldebye'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.Ldebye = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-Ldebye-cut'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.Ldebyecut = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-Ntable'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a positive integer\n')
                self.Ntable = int(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-lj'):
                if i+3 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by 4 numbers\n')
                self.eps = float(argv[i+1])
                self.sig = float(argv[i+2])
                del(argv[i:i+3])
            elif (argv[i].lower() == '-lambda'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.Lcoeff  = float(argv[i+1])
                del(argv[i:i+2])
            elif ((argv[i].lower() == '-lj-rmax'.lower()) or
                  (argv[i].lower() == '-ljrmax'.lower())):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a positive number\n')
                self.rmax_LJ = float(argv[i+1])
                del(argv[i:i+2])
            elif ((argv[i].lower() == '-lj-rmax-ratio'.lower()) or
                  (argv[i].lower() == '-ljrmaxratio'.lower())):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a positive number\n')
                self.rmax_LJ_ratio = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i] in ('-barrier-height','-barrier'
                                  '-barier-height','-barier', 'U0', 'u0')):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed ' +
                                     'by a list of numbers separated by commas (no spaces)\n')
                self.U0list = argv[i + 1].split(',')
                del(argv[i:i + 2])
            elif (argv[i].lower() in ('-charge','-charge1','-q1','-qi')):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.charge1 = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() in ('-charge2','-q2','-qj')):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.charge2 = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-rshiftLJ'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.rshiftLJ = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-rshiftQ'):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.rshiftQ = float(argv[i+1])
                del(argv[i:i+2])
            elif (argv[i].lower() == '-escale'.lower()):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a number\n')
                self.escale = float(argv[i+1])
                del(argv[i:i+2])
            elif ((argv[i].lower() == '-help') or (argv[i].lower() == '--help')
                  or
                  (argv[i].lower() == '-?') or (argv[i].lower() == '--?')):
                sys.stderr.write("\n"+usage_msg+"\n")
                exit(0)
            
            elif ((argv[i][0] == '-') and (__name__ == "__main__")):
                raise InputError('Error('+g_program_name+'):\n'+\
                    'Unrecogized command line argument "'+argv[i]+\
                    '"\n'+\
                    usage_msg)
            else:
                 i += 1

        if self.charge2 == None:
           self.charge2 = self.charge1







class GenNonBondedTables(object):
    """
    A class used for creating tabular files describing non-bonded forces
    between charged Lennard-Jones particles (with finite barrier heights).

    """

    def __init__(self):
        self.settings = GNBTSettings()

    def ParseArgs(self, argv):
        # The command above will remove arguments from argv which are
        # understood by GPSettings.ParseArgs(argv).  
        # The remaining arguments will be handled below.
        self.settings.ParseArgs(argv)


    def WriteFiles(self, outfile_lt):

        l_pair_coeffs = []
        
        # Interactions between DNA and itself:

        # Then define the interactions requested by the user:
        table_filename=self.settings.table_filename
        table_file = open(table_filename, 'w')
        calc_table_nonbonded.PrintComments(table_file)

        if self.settings.object_name != '':
            outfile_lt.write(+self.settings.object_name+' { \n'
                             '\n')

        outfile_lt.write('write_once("' +
                         self.settings.script_filename + '") {\n')

        rmax_debye = self.settings.sig + self.settings.Ldebyecut
        rmax_LJ = self.settings.rmax_LJ
        if rmax_LJ == 0.0:
            rmax_LJ = self.settings.sig*self.settings.rmax_LJ_ratio
        rmax = max(rmax_debye, rmax_LJ)
        for pair in self.settings.atom_type_pairs:
            atype1 = pair[0]
            atype2 = pair[1]
            if (atype1.find('@atom:') != 0):
                atype1 = '@atom:' + atype1
            if (atype2.find('@atom:') != 0):
                atype2 = '@atom:' + atype2
                outfile_lt.write('  pair_coeff ' + atype1 + ' ' + atype2 + 
                                 ' ' + self.settings.pair_style_name +
                                 ' ' + self.settings.table_filename +
                                 ' ' + self.settings.label+'\n')
        outfile_lt.write('}\n'
                         '\n')

        if self.settings.object_name != '':
            outfile_lt.write('\n'
                             '}  # '+self.settings.object_name+'\n'
                             '\n\n\n\n\n')

        calc_table_nonbonded.PrintTable(table_file,
                                        self.settings.eps*self.settings.escale,
                                        self.settings.sig,
                                        self.settings.Lcoeff,
                                        self.settings.rmax_LJ,
                                        self.settings.rshiftLJ,
                                        self.settings.Ldebye,
                                        self.settings.charge1,  # <-- q_i
                                        self.settings.charge2,  # <-- q_j
                                        self.settings.rshiftQ,
                                        self.settings.U0*self.settings.escale,
                                        rmax,
                                        self.settings.rmin,
                                        self.settings.Ntable,
                                        self.settings.label)

        table_file.close()




def main():
    try:
        sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str
                         +'\n\n')
        argv = [arg for arg in sys.argv]
        gen_nbt = GenNonBondedTables()
        gen_nbt.ParseArgs(argv)
        # Any remaining arguments?
        if len(argv) > 1:
            raise InputError('Error('+g_program_name+'):\n'+
                             'Unrecogized command line argument "'+argv[1]+
                             '"\n\n'+
                             usage_msg)
        # Convert all of this information to moltemplate (LT) format:
        outfile_lt = sys.stdout
        gen_nbt.WriteFiles(outfile_lt)
                    

    except (ValueError, InputError) as err:
        sys.stderr.write(str(err))
        #sys.stderr.write('\n (tell andrew to clean up his code)\n')
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()

