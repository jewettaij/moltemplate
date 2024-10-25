#!/usr/bin/env python3


from oplsaa2lt_utils import (bond_header, file_equivalences_header,
                   file_header, closing_stuff, nb_header,
                   angle_header, dihedral_header, improper_header)

from oplsaa2lt_classes import Atom, Bond, Angle, Dihedral, Improper

FILE_WITH_TYPES_AND_DIHEDRALS = "./Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.sup-2.txt"
FILE_WITH_BOND_AND_ANGLES = "./Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.sup-3.txt"

NEW_LT_FILENAME = "oplsaa2023.lt"

TYPE_CONVERSION_TUPLES = (
    ('C:', "C°"), # moltemplate doesn't like :
    ("C$", "C^"), # $ is used by moltemplate
    ("N$", "N^"), # $ is used by moltemplate
    ("O$", "O^"), # $ is used by moltemplate
    ("C#", "C|"), # moltemplate doesn't like #
    ("N*", "N§"), # leaving N* would create a mess in bonded interactions

)


def rename_type(ty: str) -> str:
    """
    A function to "clean" an atom type name, changing problematic characters

    Args:
        ty (str): the type string that has to be renamed

    Returns:
        the type name, "cleaned up"
    """
    ty = ty.strip()
    for orig, changed in TYPE_CONVERSION_TUPLES:
        ty = ty.replace(orig, changed)
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
            comment = l[DEFINITION_END+1:].strip()
            #if dihed_definition == "HC-CT-CT-C(O)"
            if "improper" in comment:
                loaded_impropers.append(
                    Improper(types=types, v1=v1, v2=v2, v3=v3, v4=v4, comment=comment))
            else:
                loaded_dihedrals.append(
                    Dihedral(types=types, v1=v1, v2=v2, v3=v3, v4=v4, comment=comment))
    return loaded_dihedrals, loaded_impropers


################################################################################################
################################################################################################
# LET'S START PARSING THE FF FILES
################################################################################################
################################################################################################

with open(FILE_WITH_BOND_AND_ANGLES) as f:
    lines = f.readlines()
bonds, angles = get_bonds_and_angles(lines)

with open(FILE_WITH_TYPES_AND_DIHEDRALS) as f:
    lines = f.readlines()
dihedrals, impropers = get_dihedrals_and_impropers(lines)

# NOTE: PROBLEM WITH SAME-TYPES BONDED INTERACTIONS...
# the same dihedral (from atom types POV) can have != parameters, based on comment...
#  e.g.:
#   102   0.000     5.500     0.00      0.0        O -C -OH-HO     carboxylic acid - aliphatic
#   210   0.000     5.000     0.00      0.0        O -C -OH-HO     benzoic acids
# the same problem is observed also for bonds and angles...
def check_uniqueness(
        of_what: list[Dihedral]|list[Improper]|list[Angle]|list[Bond],
        skip_equal_parameters: bool = False):
    for idx, it1 in enumerate(of_what):
        if it1.to_skip:
            continue
        for it2 in of_what[idx+1:]:
            if it1.types == it2.types:
                it2.to_comment = True
                it1_coeff = it1.coeff_line.lstrip("#").split("#")[0]
                it2_coeff = it2.coeff_line.lstrip("#").split("#")[0]
                #print(it1_coeff, it2_coeff)
                if it1_coeff == it2_coeff and skip_equal_parameters:
                    it2.to_skip = True

for interactions in [bonds, angles, dihedrals, impropers]:
    interactions.sort(key=lambda k: k.typename)
    interactions.sort(key=lambda k: k.sort_key)
    check_uniqueness(interactions, skip_equal_parameters=True)


################################################################################################
# LET'S ADD SOME WATER MODELS 
################################################################################################

STARTING_WAT_TYPE = 9999
wat_atoms: list[Atom] = []

# TIP3P water
# the same bonded interactions are good for TIP4P and TIP5P
# (LAMMPS proposes these k values if one want to go with a flexible TIP3P model;
#  note that they are different from the TIP3F, TIP4F and TIP5F parameters in the old oplsaa.lt)
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-0, atomic_number=16, type_str="tipO", charge="-0.830", sigma="3.188", epsilon="0.102", comment="TIP3P/F water O, long-range Coulombic solver"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-1, atomic_number=1, type_str="tipH", charge="+0.415", sigma="0.0", epsilon="0.0", comment="TIP3P/F water H, long-range Coulombic solver"))
bonds.append(Bond(types=["tipO", "tipH"], k="450.00", eq="0.9572", comment="TIP3/4/5P/F O-H"))
angles.append(Angle(types=["tipH", "tipO", "tipH"], k="55.00", eq="104.52", comment="TIP3/4/5P/F H-O-H"))

# TIP4P water
# user should change the pair_style to the one that treat internally the O-M interaction,
#   and so the O-M distance (0.1250) should be added there and not as a bond...
# also, this should not be used without fix shake, so no flexible variant
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-2, atomic_number=16, type_str="tipO", charge="0.00", sigma="3.16435", epsilon="0.16275", comment="TIP4P water O, long-range Coulombic solver"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-3, atomic_number=1, type_str="tipH", charge="+0.5242", sigma="0.0", epsilon="0.0", comment="TIP4P water H, long-range Coulombic solver"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-4, atomic_number=0, type_str="tipM", charge="-1.0484", sigma="1.0", epsilon="0.0", comment="TIP4P water M, long-range Coulombic solver"))
# bonds.append(Bond(types=["tipO", "tipM"], k="900.00", eq="0.15", comment="TIP4P O-M"))
# angles.append(Angle(types=["tipH", "tipO", "tipM"], k="50.00", eq="52.26", comment="TIP4P H-O-M"))

# TIP5P water
# user should be running this with fix rigid, so no flexible variant is provided;
#   also, bonds shouldn't matter for this model,
#   as it is kept rigid but "fix rigid" and not by bonded interactions
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-5, atomic_number=16, type_str="tipO", charge="0.00", sigma="3.0970", epsilon="0.1780", comment="TIP5P water O, long-range Coulombic solver"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-6, atomic_number=1, type_str="tipH", charge="+0.241", sigma="1.0", epsilon="0.0", comment="TIP5P water H, long-range Coulombic solver"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-7, atomic_number=0, type_str="tipL", charge="-0.241", sigma="1.0", epsilon="0.0", comment="TIP5P water L, long-range Coulombic solver"))
# bonds.append(Bond(types=["tipO", "tipL"], k="900.00", eq="0.70", comment="TIP5P O-L"))
# angles.append(Angle(types=["tipL5", "tipO", "tipL5"], k="50.00", eq="109.47", comment="TIP5P L-O-L"))
# angles.append(Angle(types=["tipH", "tipO", "tipL5"], k="50.00", eq="110.6948", comment="TIP5P H-O-L"))

# SPC and SPC/E (the same, just changes the charges on H and O...)
# should be used with fix shake, LAMMPS doesn't mention a flexible variant
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-8, atomic_number=16, type_str="spcO", charge="-0.820", sigma="3.166", epsilon="0.1553", comment="SPC water O"))
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-10, atomic_number=16, type_str="spcO", charge="-0.8476", sigma="3.166", epsilon="0.1553", comment="SPC/E water O"))
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
wat_atoms.append(Atom(type_id=STARTING_WAT_TYPE-12, atomic_number=16, type_str="opcO", charge="0.00", sigma=sigma_opc_o, epsilon="0.21280", comment="OPC water O"))
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
with open(NEW_LT_FILENAME, "w") as f:
    f.write(file_header)
    f.write("  # NOTE2: I tried to maintain the same two-letter 'general' types as from\n")
    f.write("  #  the original FF file. However, some changes had to be made to comply\n")
    f.write("  #  to the inner functioning of moltemplate. Such changes were:\n  #\n")
    for original, new in TYPE_CONVERSION_TUPLES:
        f.write(f"  #    {original} --> {new}\n")

    f.write("\n  # NOTE3: The original FF file had types for different water models,\n")
    f.write("  #  but it was missing the relevant bonded interactions; therefore, I\n")
    f.write("  #  skipped the water types from the original FF, and hardcoded some simple\n")
    f.write("  #  water models, with the relevant bonded parameters\n")

    f.write("\n  # NOTE4: Water TIP*/SPC* models parameters are taken from LAMMPS doc,\n")
    f.write("  #  the user is invited to read the proper sections in the LAMMPS user manual\n")
    f.write("  #  to properly understand how to setup a simulation with the desided model.\n")
    f.write("  #  As for OPC, it seems it could be implemented in LAMMPS similarly to the\n")
    f.write("  #   TIP4P model (where OM distance should be 0.1594 angstrom).\n")

    f.write('\n\n  write_once("In Charges") {\n')
    for line in lines[2:]:
        if line.startswith("#    Add more charge and L-J parameters"):
            break
        if line.startswith("#"):
            f.write(f"    {line}")
        else:
            atom = parse_atom_line(line)
            if atom is not None:
                atoms.append(atom)
                f.write(f"    {atom.charge_line}")
    for atom in wat_atoms:
        f.write(f"    {atom.charge_line}")
    f.write("  } # (end of atom partial charges)\n")
    f.write('\n\n  write_once("Data Masses") {\n')
    atoms += wat_atoms
    for atom in atoms:
        f.write(f"    {atom.mass_line}")
    f.write("  } # (end of atom masses)\n")

    f.write(file_equivalences_header)
    for atom in atoms:
        f.write(f"  {atom.repl_line}")

    f.write(nb_header)
    f.write('  write_once("In Settings") {\n')
    for atom in atoms:
        f.write(f"    {atom.nb_line}")
    f.write("  } # (end of pair_coeffs)\n")

    f.write("\n\n\n\n")
    f.write("  # NOTE: all bonded interaction name can't have '*' or '?' characters, so in each\n")
    f.write("  #   bonded sections such characters will be replaced with another character\n")
    f.write("  #   that, at the time of writing, is not used for atom types (* -> £, ? -> €).\n\n")
    f.write(bond_header)
    f.write('\n  write_once("In Settings") {\n')
    for bond in bonds:
        if not bond.to_skip:
            f.write(f"    {bond.coeff_line}")
    f.write("  } # (end of bond_coeffs)\n")
    f.write('\n  write_once("Data Bonds By Type") {\n')
    for bond in bonds:
        if not bond.to_skip:
            f.write(f"    {bond.bytype_line}")
    f.write("  } # (end of bonds by type)\n")

    f.write(angle_header)
    f.write('\n  write_once("In Settings") {\n')
    for angle in angles:
        if not angle.to_skip:
            f.write(f"    {angle.coeff_line}")
    f.write("  } # (end of angle_coeffs)\n")
    f.write('\n  write_once("Data Angles By Type") {\n')
    for angle in angles:
        if not angle.to_skip:
            f.write(f"    {angle.bytype_line}")
    f.write("  } # (end of angles by type)\n")

    f.write(dihedral_header)
    f.write('\n  write_once("In Settings") {\n')
    for dihedral in dihedrals:
        if not dihedral.to_skip:
            f.write(f"    {dihedral.coeff_line}")
    f.write("  } # (end of dihedral_coeffs)\n")
    f.write('\n  write_once("Data Dihedrals By Type") {\n')
    for dihedral in dihedrals:
        if not dihedral.to_skip:
            f.write(f"    {dihedral.bytype_line}")
    f.write("  } # (end of dihedrals by type)\n")

    f.write(improper_header)
    f.write('\n  write_once("In Settings") {\n')
    for improper in impropers:
        if not improper.to_skip:
            f.write(f"    {improper.coeff_line}")
    f.write("  } # (end of improper_coeffs)\n")
    f.write('\n  write_once("Data Impropers By Type (opls_imp.py)") {\n')
    for improper in impropers:
        if not improper.to_skip:
            f.write(f"    {improper.bytype_line}")
    f.write("  } # (end of impropers by type)\n")

    f.write(closing_stuff)
    f.write("}\n")
