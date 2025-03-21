#!/usr/bin/env python3

import argparse
import sys

from oplsaa2lt_utils import (bond_header, file_equivalences_header,
                   file_header, closing_stuff, nb_header,
                   angle_header, dihedral_header, improper_header)

from oplsaa2lt_classes import Atom, Bond, Angle, Dihedral, Improper

from collections import defaultdict
from distutils.version import StrictVersion

# This code only works with python 3.7 and later
version = sys.version_info
assert version.major > 3 or (version.major == 3 and version.minor >= 7)

__author__ = "Domenico Marson and Andrew Jewett"
__version__ = '1.0.1'
__date__ = '2025-2-01'

g_program_name = __file__.split('/')[-1]


TYPE_CONVERSION_TUPLES = (
    ('C:', "C°"), # moltemplate doesn't like :
    ("C$", "C^"), # $ is used by moltemplate
    ("N$", "N^"), # $ is used by moltemplate
    ("O$", "O^"), # $ is used by moltemplate
    ("C#", "C|"), # moltemplate doesn't like #
    ("N*", "N§"), # leaving N* would create a mess in bonded interactions
    ("C(O)", "C⟮"),  # moltemplate doesn't like "(" or ")", but "⟮" is okay
)


def rename_type(ty: str) -> str:
    """
    A function to "clean" an atom type name, changing problematic characters,
    or add missing characters (to ensure the type string is 2-characters long).

    Args:
        ty (str): the type string that has to be renamed

    Returns:
        the type name, "cleaned up"
    """
    ty = ty.strip()
    for orig, changed in TYPE_CONVERSION_TUPLES:
        ty = ty.replace(orig, changed)
    if ty in ["X", "Y", "Z"]:
        # In BOSS files, "X", "Y", "Z" indicate wildcards
        ty = "??"  # But here we use "??" to indicate wildcards
    if len(ty) == 1:
        # It is a good idea to make sure these strings are the same length (2).
        ty += "~" # (eg "C"->"C~")  It doesn't matter which character we add
    assert len(ty) == 2  # The "ty" strings should all be 2-characters long.
    return ty


def parse_atom_line(atom_line) -> Atom|None:
    if not atom_line[5:7].strip():
        return None
    typ, an, at, charge, sigma, epsilon, *comment = atom_line.split()
    comment = " ".join(comment)
    if "water" in comment.lower():
        return None
    return Atom(type_id=int(typ),
                atomic_number=int(an),
                type_str=rename_type(at),
                charge=charge,
                sigma=sigma,
                epsilon=epsilon,
                comment=comment)


def parse_bond_line(bond_line) -> Bond|None:
    types = list(map(rename_type, bond_line[:5].split('-')))
    k, eq, *comment = bond_line[5:].split()
    comment = " ".join(comment)
    if "water" in comment.lower() or "OW" in types or "HW" in types:
        return None
    return Bond(types=types, k=k, eq=eq, comment=comment)


def parse_angle_line(angle_line) -> Angle|None:
    types = list(map(rename_type, angle_line[:8].split('-')))
    k, eq, *comment = angle_line[8:].split()
    comment = " ".join(comment)
    if "water" in comment.lower() or "OW" in types or "HW" in types:
        return None
    return Angle(types=types, k=k, eq=eq, comment=comment)


def get_bonds_and_angles(input_lines) -> tuple[list[Bond], list[Angle]]:
    loaded_bonds: list[Bond] = []
    loaded_angles: list[Angle] = []
    read_bond = True
    for l in input_lines:
        if l.startswith("*") or l.startswith("#") or l.startswith('"""'):
            continue
        if l=="\n":
            read_bond = False
            continue
        if read_bond:
            parsed_bond = parse_bond_line(l)
            if parsed_bond is not None:
                loaded_bonds.append(parsed_bond)
        else:
            parsed_angle = parse_angle_line(l)
            if parsed_angle is not None:
                loaded_angles.append(parsed_angle)
    return loaded_bonds, loaded_angles


def get_dihedrals_and_impropers(input_lines) -> tuple[list[Dihedral], list[Improper]]:
    read_dihedrals = False
    loaded_dihedrals: list[Dihedral] = []
    loaded_impropers: list[Improper] = []
    DEFINITION_START, DEFINITION_END = 47, 60
    # NOTE: there is a strange dihedral definition; it has as the last type C(O),
    # which is not defined anywhere, so it will never be matched...
    # I kept the range 47:60, so the strange dihedral is still read (and commented), but
    # doing so creates a problem in a line in which the comments are closer to the
    # dihedral definition, hence the 'if clause' for the "P in" case.
    for l in input_lines:
        if l.startswith("#") or l.startswith("Type"):
            continue
        if "The following are the Fourier coefficients" in l:
            read_dihedrals = True
            continue
        if read_dihedrals:
            items = l[:DEFINITION_START].split()
            if len(items) <= 1:
                continue
            _type_id, v1, v2, v3, v4 = items
            dihed_definition = l[DEFINITION_START:DEFINITION_END].strip()
            if dihed_definition in ("Dummy", "", "Harmonic Rest"):
                continue
            if dihed_definition.endswith("P in"):
                dihed_definition = dihed_definition.replace("P in", "P")
            types = list(map(rename_type, dihed_definition.split('-')))
            assert len(types) == 4
            comment = l[DEFINITION_END+1:].strip()
            if "improper" in comment:
                # I've had good results using the "cenIsortJKL.py" symmetry
                # rules with OPLSAA.  (When I use those settings, the resulting
                # LAMMPS data files agrees with the data files from LigParGen.)
                # But in order for this to work, we need to swap the first two atom
                # types, since now the first atom is the center. -Andrew 2024-12-04
                types = [types[1], types[0], types[2], types[3]] # see above comment
                improper = Improper(types=types, v1=v1, v2=v2, v3=v3, v4=v4, comment=comment)
                # There is a weird problem with "allenes improper" interactions
                # that appear in the 2024 version of the BOSS files.
                # Those files specify improper interactions between atoms
                # that are typically collinear.  That would be numerically
                # unstable so we must comment these out.  This seems like
                # a bug in the BOSS files.  Perhaps later I'll report this.
                # But for now, I just comment them out.
                if comment.strip() == "allenes improper":
                    improper.to_comment = True

                loaded_impropers.append(improper)

            else:
                loaded_dihedrals.append(
                    Dihedral(types=types, v1=v1, v2=v2, v3=v3, v4=v4, comment=comment))
    return loaded_dihedrals, loaded_impropers


def delete_redundant_duplicates(
    interactions: list[Dihedral]|list[Improper]|list[Angle]|list[Bond],
) -> None:
    """ OPLSAA files contain many duplicates.  Decide which interactions to keep and return the list to the caller."""

    ##### My apologies for this ugly code. #####
    # It wasn't planned that way.  The PAR and SB files
    # that store OPLSAA parameters have grown messy over time.
    # We keep finding more redundancy problems that need fixing.
    # That's why this code has so many edge cases.

    # ----------- Step 1: ---------------
    # Organize interactions according to:
    # -atom type strings (a tuple which may include wildcards)
    # -force-field parameters  (This is called "params" in the code.)
    # -the "paramstr" (This is the entire string following the atom type list.
    #   including the parameters AND comments, if present.)
    # The result is stored in a dictionary-of-dictionary-of-dictionaries.
    # named "types_to_paraminteractions".  Lookup interactions this way:
    # types_to_paraminteractions[atomtypes][params][paramstr] --> interaction
    #
    # Later, during Step2, we will use it to figure out which interactions
    # are redundant and can be discarded.
    types_to_paraminteractions = {}
    for interaction in interactions:
        types = tuple(interaction.types)
        # Find the text containing the force field parameters ("paramstr")
        paramstr = interaction.coeff_line
        # Typically, this is "angle_coeff @angle:C3_N~_H~ 38. 118.4 # WJ94"
        # But we only want "38. 118.4 # WJ94".
        # So we ignore the text before the first two spaces:
        ispace1 = paramstr.find(" ")
        ispace2 = paramstr.find(" ", ispace1+1)
        ispace2 = ispace2+1
        paramstr = paramstr[ispace2:].strip()  # eg. "38. 118.4 # WJ94"
        # Finally, strip off the comment, leaving only the parameters
        params = paramstr.lstrip("#").split("#")[0].strip()  # eg. "38. 118.4"
        if types not in types_to_paraminteractions:
            # Multiple interactions can exist for these same atom types
            # differing by either the parameters, or comments, or both.
            # Store a dictionary that looks up the interaction between these
            # atoms according to their params (without comments).
            types_to_paraminteractions[types] = {}  # We will use this below

        # Create a dictionary that looks up all the interactions
        # which share the same types and parameters (.ie params), but may
        # have different comments.  We might want to delete these duplicates
        # later, but for now, we keep track of all of them.
        params_to_interactions = types_to_paraminteractions[types]
        if params not in params_to_interactions:
            params_to_interactions[params] = {}
        paramstr_to_interactions = params_to_interactions[params]

        # For these interactions, its possible that multiple interactions
        # may exist in the file with identical atom types, parameters (coeffs)
        # AND comments.  To get rid of these trivial duplicates, we store them in
        # a dictionary, indexed by the original paramstr (including comments).
        paramstr_to_interactions[paramstr] = interaction

    # ----------- Step 2: ---------------
    # Decide which of these interactions should be discarded.
    del interactions[:]  # we store the interactions that aren't discarded here
    for types, params_to_interactions in types_to_paraminteractions.items():
        for params, paramstr_to_interactions in params_to_interactions.items():
            for paramstr, interaction in paramstr_to_interactions.items():
                discarded = True
                paramstr_blank_comment = paramstr
                i_comment = paramstr.find("#")
                if i_comment > 0:
                    paramstr_blank_comment = paramstr[:i_comment+1]
                    comment = paramstr[i_comment+1:]
                if len(paramstr_to_interactions) == 1:
                    interactions.append(interaction)
                    discarded = False
                elif (
                    paramstr not in (params, paramstr_blank_comment)  # case 1
                    and (comment.strip() != '"')  # case 2  (see below)
                ):
                    # Edge Case 1:
                    # We want to ignore interactions if
                    # they lack a comment but they are otherwise identical.
                    # ...so we also check to make sure that paramstr != params.
                    #   (and also paramstr != paramstr_blank_comment)
                    #
                    # Edge Case 2:
                    # For some reason, a lot of comments only contain '"'.  We
                    # want to skip those too (since they are otherwise identical)
                    #
                    # If none of these edge-cases are true,
                    # then we don't discard the interaction.
                    interactions.append(interaction)
                    discarded = False
                # If all of the duplicate interactions are identical except for
                # the comments following the params, then throw away these
                # duplicates (merge them into a single interaction).
                if (len(params_to_interactions) == 1) and (not discarded):
                    break  # skip the remaining duplicates for these atoms


# NOTE: PROBLEM WITH SAME-TYPES BONDED INTERACTIONS...
# The same dihedral (from atom types POV) can have != parameters, based on comment...
#  e.g.:
#   102   0.000     5.500     0.00      0.0        O -C -OH-HO     carboxylic acid - aliphatic
#   210   0.000     5.000     0.00      0.0        O -C -OH-HO     benzoic acids
# the same problem is observed also for bonds and angles...
def count_nonredundant_duplicates(
    interactions: list[Dihedral]|list[Improper]|list[Angle]|list[Bond],
) -> None:
    for idx, it1 in enumerate(interactions):
        for it2 in interactions[idx+1:]:
            if it1.types == it2.types:
                if it1.duplicate_count == 0:
                    it1.duplicate_count = 1
                it2.duplicate_count = it1.duplicate_count + 1




def sort_duplicates(
    interactions: list[Dihedral]|list[Improper]|list[Angle]|list[Bond],
) -> None:
        # According to William Jorgensen, if duplicate bonded interactions
        # are defined for -exactly- the same atom types, then
        # the most general interaction is the one that appears earliest in
        # the PAR file (ie. the one with the lowest self.duplicate_count).
        # Since, in this case, we lack any other criteria to choose from, we
        # want moltemplate to select the most general interaction by default
        # for those atoms.  But in order to make it do that, it must appear
        # last in the generated .lt file.  (Moltemplate gives highest priority
        # to the last entry in the "By Type" list which matches the atom types.)
        # So we sort the interactions with multiple duplicates
        # according to the duplicate_count (in reverse order).
        interaction_orderkey_pairs = []
        counter = 0
        for x in interactions:
            if x.duplicate_count <= 1:
                counter += 1  # (counter is the same for each duplicate)
            k=(counter, -x.duplicate_count)  # key to be used for sorting
            interaction_orderkey_pairs.append((k, x))
        # Now sort by k
        # (ie. sort by -duplicate_count, but only when duplicate_count is > 0)
        interaction_orderkey_pairs.sort(key = lambda x: x[0])  # sort by k
        # Strip off the key and store the result in interactions.
        for i in range(len(interactions)):
            interactions[i] = interaction_orderkey_pairs[i][1]




def main(argv):
    sys.stderr.write(g_program_name + ", version " +
                     __version__ + ", " + __date__ + "\n")
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--par",
        dest="infile_name_par",
        required=True,
        # default="../oplsaa2024_original_format/Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.sup-2.par",
        help="name of the BOSS .par file.  (Example: \"oplsaa.par\")",
    )
    ap.add_argument(
        "--sb",
        dest="infile_name_sb",
        required=True,
        # default="../oplsaa2024_original_format/Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.sup-3.sb",
        help="name of the BOSS .sb file.  (Example: \"oplsaa.sb\")",
    )
    ap.add_argument(
        "--out",
        dest="outfile_name",
        required=True,
        help="name of the .lt file you want to create.",
    )
    ap.add_argument(
        "--name",
        dest="object_name",
        default="OPLSAA",
        required=False,
        help='name of the moltemplate object you want to create. (Default: "OPLSAA")',
    )

    args = ap.parse_args(argv[1:])

    if args.outfile_name != "":
        outfile = open(args.outfile_name, 'w')
    else:
        outfile = sys.stdout  # print to the terminal

    ################################################################################################
    ################################################################################################
    # LET'S START PARSING THE FF FILES
    ################################################################################################
    ################################################################################################

    # Import bond and angle information from the ".sb" BOSS file.
    # (Eg. https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.3c06602/suppl_file/jp3c06602_si_003.txt)
    with open(args.infile_name_sb, "r") as infile_sb:
       lines_sb = infile_sb.readlines()
    bonds, angles = get_bonds_and_angles(lines_sb)


    # Import atom, dihedral, and improper information from the ".par" BOSS file.
    # (Eg. https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.3c06602/suppl_file/jp3c06602_si_002.txt)
    with open(args.infile_name_par, "r") as infile_par:
       lines_par = infile_par.readlines()
    dihedrals, impropers = get_dihedrals_and_impropers(lines_par)


    # Now, let's cleanup all the lists of bonded interactions.
    # (Remove redundant entries, and sort by atom-type)
    for interactions in [bonds, angles, dihedrals, impropers]:
        interactions.sort(key=lambda x: x.typename)
        interactions.sort(key=lambda x: x.sort_key)
        delete_redundant_duplicates(interactions)
        count_nonredundant_duplicates(interactions)
        sort_duplicates(interactions)


    ################################################################################################
    # LET'S ADD SOME WATER MODELS 
    ################################################################################################

    STARTING_WAT_TYPE = 9999
    wat_atoms: list[Atom] = []

    # TIP3P water
    # Note: TIP3P shares the same bonded interactions with TIP4P, TIP5P, etc...
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-0, atomic_number=8, type_str="tipO", charge="-0.834", sigma="3.188", epsilon="0.102", comment="TIP3P water O, long-range Coulombic solver"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-1, atomic_number=1, type_str="tipH", charge="+0.417", sigma="0.0", epsilon="0.0", comment="TIP3P water H, long-range Coulombic solver"))
    # (LAMMPS proposes using these k values if you want to use a flexible TIP3P model.)
    bonds.append(Bond(types=["tipO", "tipH"], k="450.00", eq="0.9572", comment="TIP3/4/5P/F O-H"))
    angles.append(Angle(types=["tipH", "tipO", "tipH"], k="55.00", eq="104.52", comment="TIP3/4/5P/F H-O-H"))

    #############################################
    # COMMENTING OUT: I can't figure out which water model this refers to: TIP3P/F
    # TIP3P/F water   <--Is it TIP3P/Fs?  (If so, the charges are wrong.)
    # wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-0, atomic_number=8, type_str="tipO", charge="-0.830", sigma="3.188", epsilon="0.102", comment="TIP3P/F water O, long-range Coulombic solver"))
    # wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-1, atomic_number=1, type_str="tipH", charge="+0.415", sigma="0.0", epsilon="0.0", comment="TIP3P/F water H, long-range Coulombic solver"))
    # bonds.append(Bond(types=["tipO", "tipH"], k="450.00", eq="0.9572", comment="TIP3/4/5P/F O-H"))
    # angles.append(Angle(types=["tipH", "tipO", "tipH"], k="55.00", eq="104.52", comment="TIP3/4/5P/F H-O-H"))
    #############################################

    # TIP4P water
    # user should change the pair_style to the one that treat internally the O-M interaction,
    #   and so the O-M distance (0.1250) should be added there and not as a bond...
    # also, this should not be used without fix shake, so no flexible variant
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-2, atomic_number=8, type_str="tipO", charge="0.00", sigma="3.16435", epsilon="0.16275", comment="TIP4P water O, long-range Coulombic solver"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-3, atomic_number=1, type_str="tipH", charge="+0.5242", sigma="0.0", epsilon="0.0", comment="TIP4P water H, long-range Coulombic solver"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-4, atomic_number=0, type_str="tipM", charge="-1.0484", sigma="1.0", epsilon="0.0", comment="TIP4P water M, long-range Coulombic solver"))
    # bonds.append(Bond(types=["tipO", "tipM"], k="900.00", eq="0.15", comment="TIP4P O-M"))
    # angles.append(Angle(types=["tipH", "tipO", "tipM"], k="50.00", eq="52.26", comment="TIP4P H-O-M"))

    # TIP5P water
    # user should be running this with fix rigid, so no flexible variant is provided;
    #   also, bonds shouldn't matter for this model,
    #   as it is kept rigid but "fix rigid" and not by bonded interactions
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-5, atomic_number=8, type_str="tipO", charge="0.00", sigma="3.0970", epsilon="0.1780", comment="TIP5P water O, long-range Coulombic solver"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-6, atomic_number=1, type_str="tipH", charge="+0.241", sigma="1.0", epsilon="0.0", comment="TIP5P water H, long-range Coulombic solver"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-7, atomic_number=0, type_str="tipL", charge="-0.241", sigma="1.0", epsilon="0.0", comment="TIP5P water L, long-range Coulombic solver"))
    # bonds.append(Bond(types=["tipO", "tipL"], k="900.00", eq="0.70", comment="TIP5P O-L"))
    # angles.append(Angle(types=["tipL5", "tipO", "tipL5"], k="50.00", eq="109.47", comment="TIP5P L-O-L"))
    # angles.append(Angle(types=["tipH", "tipO", "tipL5"], k="50.00", eq="110.6948", comment="TIP5P H-O-L"))

    # SPC and SPC/E (the same, just changes the charges on H and O...)
    # should be used with fix shake, LAMMPS doesn't mention a flexible variant
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-8, atomic_number=8, type_str="spcO", charge="-0.820", sigma="3.166", epsilon="0.1553", comment="SPC water O"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-10, atomic_number=8, type_str="spcO", charge="-0.8476", sigma="3.166", epsilon="0.1553", comment="SPC/E water O"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-9, atomic_number=1, type_str="spcH", charge="+0.410", sigma="0.0", epsilon="0.0", comment="SPC water H"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-11, atomic_number=1, type_str="spcH", charge="+0.4238", sigma="0.0", epsilon="0.0", comment="SPC/E water H"))
    bonds.append(Bond(types=["spcO", "spcH"], k="450.00", eq="1.000", comment="SPC-SPC/E O-H"))
    angles.append(Angle(types=["spcH", "spcO", "spcH"], k="55.00", eq="109.47", comment="SPC-SPC/E H-O-H"))

    # OPC
    # I saw users using this water model via the lj/long/tip4p/long pair_style,
    #   so as for TIP4P users need to provide O-E distance (0.1594) there and use fix shake.
    # parameters are taken from AmberTools2024
    #   where LJ distance parameter is provided as half r_min, not sigma, hence the conversion
    HALFRMIN2SIGMA = 2/(2**(1/6))
    sigma_opc_o = f"{1.777167268 * HALFRMIN2SIGMA:10.6f}".strip()
    sigma_opc_ep = f"{1 * HALFRMIN2SIGMA:10.6f}".strip()
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-12, atomic_number=8, type_str="opcO", charge="0.00", sigma=sigma_opc_o, epsilon="0.21280", comment="OPC water O"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-13, atomic_number=1, type_str="opcH", charge="+0.679142", sigma="0.0", epsilon="0.0", comment="OPC water H"))
    wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-14, atomic_number=0, type_str="opcE", charge="-1.358284", sigma=sigma_opc_ep, epsilon="0.0", comment="OPC water E"))
    bonds.append(Bond(types=["opcO", "opcH"], k="450.00", eq="0.8724", comment="OPC O-H"))
    angles.append(Angle(types=["opcH", "opcO", "opcH"], k="55.00", eq="103.6", comment="OPC H-O-H"))
    # bonds.append(Bond(types=["opcO", "opcE"], k="450.00", eq="0.1594", comment="OPC O-LP"))



    ################################################################################################
    ###### LET'S START WRITING STUFF
    ##################### NOTE: i load and write the atom types concurrently,
    #####################       in this way I can keep the comment blocks from the original FF file
    ################################################################################################
    atoms: list[Atom] = []
    with open(args.outfile_name, "w") as outfile:
        outfile.write("# This file was generated automatically using:\n")
        outfile.write("# " + g_program_name + " " + " ".join(argv[1:]) + "\n")
        outfile.write("")
        outfile.write(file_header)
        outfile.write("  # NOTE2: I tried to maintain the same two-letter 'general' types as from\n")
        outfile.write("  #  the original FF file. However, some changes had to be made to comply\n")
        outfile.write("  #  to the inner functioning of moltemplate. Such changes were:\n  #\n")
        for original, new in TYPE_CONVERSION_TUPLES:
            outfile.write(f"  #    {original} --> {new}\n")

        outfile.write("\n  # NOTE3: The original FF file had types for different water models,\n")
        outfile.write("  #  but it was missing the relevant bonded interactions; therefore, I\n")
        outfile.write("  #  skipped the water types from the original FF, and hardcoded some simple\n")
        outfile.write("  #  water models, with the relevant bonded parameters\n")

        outfile.write("\n  # NOTE4: Water TIP*/SPC* models parameters are taken from LAMMPS doc,\n")
        outfile.write("  #  the user is invited to read the proper sections in the LAMMPS user manual\n")
        outfile.write("  #  to properly understand how to setup a simulation with the desided model.\n")
        outfile.write("  #  As for OPC, it seems it could be implemented in LAMMPS similarly to the\n")
        outfile.write("  #   TIP4P model (where OM distance should be 0.1594 angstrom).\n")

        outfile.write('\n\n  write_once("In Charges") {\n')
        for line in lines_par[2:]:
            if line.startswith("#    Add more charge and L-J parameters"):
                break
            if line.startswith("#"):
                outfile.write(f"    {line}")
            else:
                atom = parse_atom_line(line)
                if atom is not None:
                    atoms.append(atom)
                    outfile.write(f"    {atom.charge_line}")
        for atom in wat_atoms:
            outfile.write(f"    {atom.charge_line}")
        outfile.write("  } # (end of atom partial charges)\n")
        outfile.write('\n\n  write_once("Data Masses") {\n')
        atoms += wat_atoms
        for atom in atoms:
            outfile.write(f"    {atom.mass_line}")
        outfile.write("  } # (end of atom masses)\n")

        outfile.write(file_equivalences_header)
        for atom in atoms:
            outfile.write(f"  {atom.repl_line}")

        outfile.write(nb_header)
        outfile.write('  write_once("In Settings") {\n')
        for atom in atoms:
            outfile.write(f"    {atom.nb_line}")
        outfile.write("  } # (end of pair_coeffs)\n")

        outfile.write("\n\n\n\n")
        outfile.write("  # NOTE: all bonded interaction name can't have '*' or '?' characters, so in each\n")
        outfile.write("  #   bonded sections such characters will be replaced with another character\n")
        outfile.write("  #   that, at the time of writing, is not used for atom types (* -> £, ? -> €).\n\n")
        outfile.write(bond_header)
        outfile.write('\n  write_once("In Settings") {\n')
        for bond in bonds:
            if not bond.to_skip:
                outfile.write(f"    {bond.coeff_line}")
        outfile.write("  } # (end of bond_coeffs)\n")
        outfile.write('\n  write_once("Data Bonds By Type") {\n')
        for bond in bonds:
            if not bond.to_skip:
                outfile.write(f"    {bond.bytype_line}")
        outfile.write("  } # (end of bonds by type)\n")

        outfile.write(angle_header)
        outfile.write('\n  write_once("In Settings") {\n')
        for angle in angles:
            if not angle.to_skip:
                outfile.write(f"    {angle.coeff_line}")
        outfile.write("  } # (end of angle_coeffs)\n")
        outfile.write('\n  write_once("Data Angles By Type") {\n')
        for angle in angles:
            if not angle.to_skip:
                outfile.write(f"    {angle.bytype_line}")
        outfile.write("  } # (end of angles by type)\n")

        outfile.write(dihedral_header)
        outfile.write('\n  write_once("In Settings") {\n')
        for dihedral in dihedrals:
            if not dihedral.to_skip:
                outfile.write(f"    {dihedral.coeff_line}")
        outfile.write("  } # (end of dihedral_coeffs)\n")
        outfile.write('\n  write_once("Data Dihedrals By Type") {\n')
        for dihedral in dihedrals:
            if not dihedral.to_skip:
                outfile.write(f"    {dihedral.bytype_line}")
        outfile.write("  } # (end of dihedrals by type)\n")

        outfile.write(improper_header)
        outfile.write('\n  write_once("In Settings") {\n')
        for improper in impropers:
            if not improper.to_skip:
                outfile.write(f"    {improper.coeff_line}")
        outfile.write("  } # (end of improper_coeffs)\n")
        outfile.write('\n  write_once("Data Impropers By Type (cenIsortJKL.py)") {\n')
        for improper in impropers:
            if not improper.to_skip:
                outfile.write(f"    {improper.bytype_line}")
        outfile.write("  } # (end of impropers by type)\n")

        outfile.write(closing_stuff)
        outfile.write("}\n")


if __name__ == '__main__':
    main(sys.argv)
