#!/usr/bin/env python3

# Author: Andrew Jewett
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2022


import sys
import argparse
from operator import itemgetter
from collections import defaultdict


# Global variables
g_filename = __file__.split('/')[-1]
g_module_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str = '2022-8-21'
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





def ReadCharges(fcharges):
    """ An optional function which reads charges from an external file.
        (instead of the MOL2 file).
        That file can have one of two possible formats.  1-Column format:
        -0.18
         0.06
         0.06
         0.06
          :
    """
    charges = []
    for line in fcharges:
        # First remove any commented text following a # character
        ic = line.find('#')
        if ic != -1:
            line = line[0:ic]
        # Then discard empty lines
        line = line.strip()
        if line == "":
            continue
        tokens = line.split()
        # Does the file contain one column or two columns?
        charge = float(tokens[0])
        charges.append(charge)
    return charges





def ConvertMol22Lt(fin = sys.stdin,
                   fout = sys.stdout,
                   ff_name = None,
                   ff_filename = None,
                   object_name = None,
                   fcharges = None,
                   atom_type_capitalization = None,
                   atom_name_capitalization = None,
                   include_bond_types = False):
    """ Reads a MOL2 file (and an optional file containing charges)
        and creates a moltemplate file (LT file). """
    assert(fin)
    assert(fout)

    # If each subunit is part of a larger molecule (named object_name) then
    # merge the molecule-IDs together.  All atoms are part of the same molcule.
    merge_mol_ids = (object_name != None) and (object_name != "")

    external_charges = []
    # If there is an additional file containing charges, read that first
    if (fcharges):
        external_charges = ReadCharges(fcharges)

    id2name = []     # lookup atom name from ID number
    id2type = []     # lookup atom name from ID number
    id2subId = []
    id2charge = []
    id2coords = []
    bonds_orig = []  # the list of bonds in the MOL2 file in the original order
    subId2atomNames   = defaultdict(dict)
    subId2subNameOrig = {}
    subId2subName     = {}
    subNameOrig2subId = defaultdict(list)
    subNamesUsed      = set([])

    section_name = None

    ############################################
    ############  Read the MOL2 file ###########
    ############################################

    # Now read the MOL2 file, one line at a time
    for line in fin:
        # First remove any commented text following a # character
        ic = line.find('#')
        if ic != -1:
            line = line[0:ic]
        # Then discard empty lines
        line = line.strip()
        if line == "":
            continue
        if line.find('@<TRIPOS>') == 0:
            section_name=line[9:] #=everything after "@<TRIPOS>" (eg "ATOM")
            continue
        tokens = line.split() # split line into columns
        # What we do next, depends on which section of the file we are in:

        # Read the lines in the "ATOM" section
        if section_name == "ATOM":
            if len(tokens) < 6:
                raise InputError("Error: The ATOM section of a MOL2 file should contain at least 6 columns.")
            # Parse the entries on each line
            atom_id    = int(tokens[0])
            atom_name  = tokens[1]
            x          = float(tokens[2])
            y          = float(tokens[3])
            z          = float(tokens[4])
            atom_type_name = tokens[5]
            if len(tokens) > 6:
                sub_id     = int(tokens[6])
            if len(tokens) > 7:
                sub_name   = tokens[7]
            charge     = 0
            if len(tokens) > 8:
                charge = float(tokens[8])
            if len(external_charges) > 0:
                if atom_id-1 < len(external_charges):
                    charge = external_charges[atom_id-1]
                else:
                    raise InputError("Error: The MOL2 file is longer than the charges file.")

            # Process the information on this line
            # First make sure all arrays are large enough to store the new info
            if atom_name_capitalization == 'upper':
                atom_name = atom_name.upper()
            elif atom_name_capitalization == 'lower':
                atom_name = atom_name.lower()
            if atom_type_capitalization == 'upper':
                atom_type_name = atom_type_name.upper()
            elif atom_type_capitalization == 'lower':
                atom_type_name = atom_type_name.lower()
            if atom_id >= len(id2name):
                num_ids = len(id2name)
                assert(num_ids == len(id2type))
                assert(num_ids == len(id2subId))
                id2name    += (1 + atom_id - num_ids) * [""]
                id2type    += (1 + atom_id - num_ids) * [""]
                id2subId   += (1 + atom_id - num_ids) * [-1]
                id2charge  += (1 + atom_id - num_ids) * [0]
                id2coords  += (1 + atom_id - num_ids) * [[-1.0,-1.0,-1.0]]
            # Then, make sure this atomID number is unique
            if id2name[atom_id] != "":
                raise InputError("Format Error: Repeated atom-ID number in ATOM section of MOL2 file: "+str(atom_id))
            # Now store the result
            id2name[atom_id]   = atom_name
            id2type[atom_id]   = atom_type_name
            id2subId[atom_id]  = sub_id
            id2charge[atom_id] = charge
            id2coords[atom_id][0] = x
            id2coords[atom_id][1] = y
            id2coords[atom_id][2] = z
            # Make sure the array with the info for each subunit is large enough
            subId2subNameOrig[sub_id] = sub_name
            # For convenience, rename subId2atomNames[sub_id] to "atomName2id"
            atomName2id = subId2atomNames[sub_id]
            # Atom names are not unique.
            # (For example, multiple monomers may have an atom named "C1".)
            # But aton names should be unique within each subunit.
            # Check that this is so.
            if atom_name in atomName2id:
                raise InputError('Format Error: Repeated atom-name "'+atom_name+' in subunit #'+str(sub_id)+' of MOL2 file.')
            # If this atom_name has not been encountered yet in this subunit
            # then add it to the subunit, and store the corresponding atom_id
            atomName2id[atom_name] = atom_id

        # Read the lines in the "BOND" section
        elif section_name == "BOND":
            if len(tokens) < 3:
                raise InputError("Error: The BOND section of a MOL2 file should contain at least 3 columns.")
            # Parse the entries on each line
            #bond_id     = int(tokens[0])   #discard
            atom1_id    = int(tokens[1])
            atom2_id    = int(tokens[2])
            bond_type   = ""
            if len(tokens) > 3:
                bond_type = tokens[3]
            bonds_orig.append([atom1_id, atom2_id, bond_type])

    ############################################
    ##############  Postprocessing #############
    ############################################

    # First, check for missing or non consecutive atomId numbers
    for i in range(1, len(id2name)):
        if id2name[i] == "":
            raise InputError("Format Error: Missing atom-ID in ATOM section of MOL2 file: "+str(i))
        assert(id2type[i] != "")

    # How many molecular subunits are present in this MOL2 file?
    num_subunits = len(subId2atomNames)
    assert(num_subunits == len(subId2subNameOrig))
    global_bonds  = []
    subId2subName = {}

    # If the the MOL2 file contains multiple identical types of molecules
    # or molecular subunits, the resulting LT file will contain multiple
    # redundant definitions of the same molecular subunits.
    # This won't cause any problems (other than large files).
    # If the user wants to avoid redefining the same types of molecules
    # (or molecular subunits), it is their responsibility
    # to supply a MOL2 file which does not contain this redudancy.
    # (Then they can use moltemplate's "new", ".move()", and ".rot()" commands
    # to instantiate multiple copies of the molecule instead of redefining it.)
    #
    # As mentioned above, for this program, each type of molecular subunit in
    # the MOL2 file is assumed to be unique.  In MOL2 files, multiple
    # molecular subunits can have the same name (but different subunit IDs).
    # However in the moltemplate file that we create, each type of molecular
    # subunit must have a unique name.
    # So we must loop over all the subunits and change their names to
    # make sure that every subunit has a unique name.
    # (The new name will be a combination of the original subunit name
    # and the subunit-ID number.)
    # First count how many times the same name is used for different subIds
    for sub_id, sub_name in subId2subNameOrig.items():
        assert(sub_name != "")
        subNameOrig2subId[sub_name].append(sub_id)

    # Then, for each subunit with a duplicated name, add a numeric suffix
    for sub_name, sub_ids in sorted(subNameOrig2subId.items(),
                                    key = itemgetter(1)):
        if len(sub_ids) == 1:
            sub_id = sub_ids[0]
            subId2subName[sub_id] = sub_name
            # If the caller requested a custom name for the molecule, use that.
            if (object_name and (object_name != "") and (num_subunits == 1)):
                subId2subName[sub_id] = object_name
            subNamesUsed.add(sub_name)
        else:
            assert(len(sub_ids) > 1)
            for sub_id in sub_ids:
                # It's possible that when we add a new numeric suffix, it
                # might clash with one of the names that was already used.
                # So try appending different suffixes to the end of the name
                # until we avoid a name clash.  Do this by adding a variable
                # number of '_' characters to the end of the subunit name
                # followed by the subunit number until the new name
                # hasn't been encountered already.
                n_blanks = 1
                while True:
                    new_sub_name = sub_name + '_'*n_blanks + str(sub_id)
                    if not (new_sub_name in subNamesUsed):
                        subId2subName[sub_id] = new_sub_name
                        subNamesUsed.add(new_sub_name)
                        break
                    n_blanks += 1

    ##### Now go through the list of bonds.
    #
    # Most bonds are between atoms in the same molecule or subunit.
    # But some bonds connect different molecular subunits together.
    # So we classify each bond as one of these two types:
    # Bonds between atoms in the same subunit should be grouped together
    # and written out as part of that molecular subunit's definition.
    # (Bonds between atoms in different subunits should be written out later.)
    subId2bonds   = defaultdict(list)
    for ib in range(0, len(bonds_orig)):
        ia1, ia2, bond_type = bonds_orig[ib]
        if id2subId[ia1] == id2subId[ia2]:
            subId2bonds[id2subId[ia1]].append([ia1, ia2, bond_type])
        else:
            global_bonds.append([ia1, ia2, bond_type])


    ############################################
    ### Write the moltemplate file (LT file) ###
    ############################################

    if ff_filename and ff_filename != "":
        fout.write('import "' + ff_filename + '"\n\n\n')

    ff_str = ""
    if ff_name and ff_name != "":
        ff_str = ' inherits ' + ff_name

    if (object_name and (object_name != "") and (num_subunits > 1)):
        fout.write(object_name + ff_str + ' {\n\n')
        fout.write('# (Optional: Let moltemplate know that the atoms in each molecular subunit\n'
                   '#  share the same molecule-ID "..." using the "create_var" command)\n\n')
        fout.write('create_var {$mol}\n\n')

    # Now loop over all of the types of molecular subunits
    # defined in the MOL2 file, and convert them to moltemplate-style
    # molecule definitions (each containing a "Data Atoms" and "Bonds" section)
    for sub_id, sub_name in sorted(subId2subName.items(),
                                   key = itemgetter(1)):
        # Print out the "Data Atoms" section of the LT file
        mol_id_name = "m"   # default molecule ID name
        # If the are multiple subunits which are part of the same molecule then
        # choose a global mol_id_name "..." shared by all atoms in this molecule
        if (num_subunits > 1) and merge_mol_ids:
            mol_id_name = "..." #this is moltemplate syntax for sharing counters

        ###### Print the name of this type of molecular subunit followed by {}.
        # This will create a block of text enclosed in {} parenthesis
        # which defines the type of molecular subunit.
        fout.write(sub_name + ff_str + ' {\n'
                   '\n')

        ###### Print the "Data Atoms" section for this molecular subunit #####
        fout.write('  #  atomId molId atomType charge X Y Z\n'
                   '\n'
                   '  write("Data Atoms"){\n')
        # Now loop through the atoms in each subunit, sortedy by atom_id
        for atom_name, atom_id in sorted(subId2atomNames[sub_id].items(),
                                         key = itemgetter(1)):
            fout.write('    $atom:'+str(atom_name) +
                       ' $mol:' + mol_id_name +
                       ' @atom:' + id2type[atom_id] + ' ' +
                       str(id2charge[atom_id]) + ' ' +
                       str(id2coords[atom_id][0]) + ' ' +
                       str(id2coords[atom_id][1]) + ' ' +
                       str(id2coords[atom_id][2]) + '\n')
        fout.write('  } # Atoms section\n')

        ######## Print the "Bonds" section for this molecular subunit ########
        #
        # Loop through the bonds in each subunit:
        if include_bond_types:
            fout.write('\n'
                       '  #  bondId  bondType  atomId1 atomId2\n'
                       '\n')
            fout.write('  write("Data Bonds"){\n')
        else:
            fout.write('\n'
                       '  #  bondId  atomId1 atomId2\n'
                       '\n')
            fout.write('  write("Data Bond List") {\n')
        bonds = subId2bonds[sub_id]
        for ib in range(0, len(bonds)):
            assert(len(bonds[ib]) == 3)
            ia1, ia2, btype = bonds[ib]
            fout.write('    $bond:b'+str(ib+1))
            if include_bond_types:
                fout.write('  @bond:'+str(bond_type))
            fout.write('  $atom:' + id2name[ia1] +
                       '  $atom:' + id2name[ia2])
            #if (not include_bond_types) and (btype != ""):
            #    fout.write('  # suggested bond type: ' + bond_type)
            fout.write('\n')
        fout.write('  } # Bonds section\n\n')

        # We are done with this molecular subunit's definition
        fout.write('}  # '+sub_name+'\n\n\n\n')

    # If there are multiple molecular subunits
    if ((object_name and object_name != "" and num_subunits>1) or
        (len(global_bonds) > 0)):

        fout.write('# Now instantiate a copy of each molecular subunit we defined earlier.\n\n')

        assert((len(global_bonds) > 0) == (num_subunits > 1))
        for sub_id, sub_name in sorted(subId2subName.items(),
                                       key = itemgetter(1)):
            fout.write(sub_name + '_instance = new ' + sub_name + '\n')

        if not (object_name and object_name != ""):
            usage_instructions = \
                '# -------- INSTRUCTIONS FOR USING THIS FILE: --------\n'+\
                '# You can either run moltemplate.sh directly on this file\n'+\
                '# or create a new LT file (eg "system.lt") and use moltemplate\'s\n'+\
                '# "import" command to import this command in your system.lt file.\n'+\
                '#\n'+\
                '# Alternatively, you want to make multiple copies of the atoms in this file\n'+\
                '# then re-run '+g_program_name+' with the "--name MOL_NAME" argument.\n'+\
                '# This will encapsulate all of the text in this file within a molecule object\n'+\
                '# (named MOL_NAME), which you can easily make multiple copies of later using\n'+\
                '# moltemplate\'s "new" command.  For example:\n'+\
                '# copy1 = new MOL_NAME\n'+\
                '# copy2 = new MOL_NAME.move(5.2,0,0).rot(180,1,0,0)\n'+\
                '# copy3 = new MOL_NAME.move(10.4,0,0)\n'+\
                '# copy4 = new MOL_NAME.move(15.6,0,0).rot(180,1,0,0)\n'+\
                '#   :   =   :\n'+\
                '# ---------------------------------------------------\n'

    else:
        usage_instructions = \
            '# -------- INSTRUCTIONS FOR USING THIS FILE: --------\n'+\
            '# So far, we have just defined (one or more) molecular subunits.\n' +\
            '# If you want to use these molecule(s) in a simulation, you must instantiate\n' +\
            '# copies of them.  To do that you would the "new" command.  For example:\n' +\
            '#\n'
        for sub_id, sub_name in sorted(subId2subName.items(),
                                       key = itemgetter(1)):
            usage_instructions += '# '+sub_name+'_instance' + \
                ' = new ' + sub_name + '\n'
        usage_instructions += \
            '#\n'+ \
            '# You could either put this command here, or in a separate file.\n' +\
            '# (...Such as "system.lt".  In that case remember to use moltemplate\'s\n' +\
            '#  "import" command to import this file beforehand because you must ensure\n' +\
            '#  that the molecules in this file are loaded before they are used.)\n' +\
            '#\n' +\
            '# Note: You can also modify the position and orientation of each copy\n' +\
            '# using the .move() and .rot() commands.  (See the moltemplate manual.)\n' +\
            '# ---------------------------------------------------\n'


    fout.write('\n\n')


    ####### We are done printing out the definitions of each molecular subunit
    ####### Now print out global information (not specific to any one subunt).

    # Print out the bonds that connect atoms in
    # different molecular subunits (if any)

    if len(global_bonds) > 0:
        fout.write('# Bonds between atoms in different molecular subunits\n'
                   '\n')
        if include_bond_types:
            fout.write('write("Data Bonds"){\n')
        else:
            fout.write('write("Data Bond List") {\n')
        bonds = global_bonds
        for ib in range(0, len(bonds)):
            assert(len(bonds[ib]) == 3)
            ia1, ia2, btype = bonds[ib]
            subname1 = subId2subName[id2subId[ia1]]
            subname2 = subId2subName[id2subId[ia2]]
            fout.write('  $bond:b'+str(ib+1))
            if include_bond_types:
                fout.write('  @bond:'+str(bond_type))
            fout.write('  $atom:' + subname1 + '_instance/' + id2name[ia1] +
                       '  $atom:' + subname2 + '_instance/' + id2name[ia2])
            #if (not include_bond_types) and (btype != ""):
            #    fout.write('  # suggested bond type: ' + bond_type)
            fout.write('\n')
        fout.write('} # global bonds section\n\n')


    if object_name and object_name != "" and num_subunits>1:
        fout.write('} # ' + object_name + '\n\n\n')
        usage_instructions = \
            '# -------- INSTRUCTIONS FOR USING THIS FILE: --------\n'+\
            '# So far, we have just defined a molecule named "'+object_name+'"\n'+ \
            '# If you want to use this molecule in a simulation, you must instantiate\n'+\
            '# a copy of it.  To do that you would the "new" command.  For example:\n'+\
            '#\n'+\
            '# ' + object_name + '_instance = new ' +object_name+ '\n'+\
            '#\n'+\
            '# You could either put this command here, or in a separate file.\n'+\
            '# (...Such as "system.lt".  In that case remember to use moltemplate\'s\n'+\
            '#  "import" command to import this file beforehand because you must ensure\n'+\
            '#  that the molecule in this file is loaded before it is used.)\n'+\
            '# ---------------------------------------------------\n'

    # Finally print out a comment with a suggestion how to use this file.
    fout.write('\n' + usage_instructions)

    # A final check to make sure that if the user supplied a file
    # containing custom charges, the length of that file should equal
    # the number of atoms in the MOL2 file
    if ((len(external_charges) > 0) and
        (len(external_charges) != len(id2type)-1)):
        raise InputError("Error: The MOL2 file is shorter than the charges file.")




def main():
    # Inform the user what version of the software they are using
    sys.stderr.write(g_program_name + ' v' +
                     g_version_str + ' ' + g_date_str + '\n')
    sys.stderr.write('WARNING: THIS IS EXPERIMENTAL SOFTWARE (2022-8-15)\n')
    try:
        # Now parse the arguments passed to the program (if any)
        ap = argparse.ArgumentParser()
        ap.add_argument('-i', '-in', '--in',
                        dest='mol2_filename',
                        required=False,
                        help='name of the mol2 file you want to convert (if unspecified, stdin is used)')
        ap.add_argument('-o', '-out', '--out',
                        dest='lt_filename',
                        required=False,
                        help='name of the LT file you want to create (if unspecified, stdout is used)')
        ap.add_argument('-q', '-charge', '-charges', '--charge', '--charges',
                        dest='charge_filename',
                        required=False,
                        help='name of a text file containing the charge of each atom (optional)')
        ap.add_argument('-name', '--name',
                        dest='object_name',
                        required=False,
                        help='name of the molecule or molecular subunit you want to create.')
        ap.add_argument('-ff', '--ff',
                        dest='ff_name',
                        required=False,
                        help='name of the force field you are using (eg "GAFF2")')
        ap.add_argument('-ff-file', '--ff-file',
                        dest='ff_filename',
                        required=False,
                        help='name of the file containing force field information (eg "gaff2.lt")')
        ap.add_argument('-upper-case-types', '--upper-case-types',
                        dest='upper_case_types',
                        required=False,
                        default=False,
                        action='store_true',
                        help='convert the atom type names to upper-case characters?')
        ap.add_argument('-lower-case-types', '--lower-case-types',
                        dest='lower_case_types',
                        required=False,
                        default=False,
                        action='store_true',
                        help='convert the atom type names to lower-case characters?')
        ap.add_argument('-upper-case-names', '--upper-case-names',
                        dest='upper_case_names',
                        required=False,
                        default=False,
                        action='store_true',
                        help='convert the atom type names to upper-case characters?')
        ap.add_argument('-lower-case-names', '--lower-case-names',
                        dest='lower_case_names',
                        required=False,
                        default=False,
                        action='store_true',
                        help='convert the atom type names to lower-case characters?')
        args = ap.parse_args()
        # Now figure out the names of the file(s) the user wants us to read
        # (By default, this program will read from the terminal (sys.stdin).)
        if args.mol2_filename:
            fmol2 = open(args.mol2_filename, 'r')
        else:
            fmol2 = sys.stdin
        fcharges = None
        if args.charge_filename:
            fcharges = open(args.charge_filename, 'r')
        # Now figure out the names of the LT file the user wants to create
        # (By default, this program will write to the terminal (sys.stdout).)
        if args.lt_filename:
            flt = open(args.lt_filename, 'w')
        else:
            flt = sys.stdout
        atom_type_capitalization = ''
        if args.upper_case_types:
            atom_type_capitalization = 'upper'
        elif args.lower_case_types:
            atom_type_capitalization = 'lower'
        atom_name_capitalization = ''
        if args.upper_case_names:
            atom_name_capitalization = 'upper'
        elif args.lower_case_names:
            atom_name_capitalization = 'lower'
        # Now convert the MOL2 file into an LT file
        ConvertMol22Lt(fmol2, flt,
                       ff_name = args.ff_name,
                       ff_filename = args.ff_filename,
                       object_name = args.object_name,
                       atom_type_capitalization = atom_type_capitalization,
                       atom_name_capitalization = atom_name_capitalization,
                       fcharges = fcharges)

    except (ValueError, InputError) as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(-1)


if __name__ == '__main__':
    main()
