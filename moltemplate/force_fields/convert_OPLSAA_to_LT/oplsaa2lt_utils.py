#!/usr/bin/env python3

file_header = """
# This file contains OPLSAA parameters and rules for creating angle, dihedral,
# and improper interactions according to OPLSAA conventions, taken from
# "Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.sup-1"
#
# USAGE: You can create molecules using this force-field this way:
#
# MyMolecule inherits OPLSAA {
#  # atom-id  mol-id atom-type charge    X        Y        Z
#  write('Data Atoms') {
#    $atom:C1  $mol @atom:88  0.00  -0.6695    0.000000  0.000000
#    $atom:H11 $mol @atom:89  0.00  -1.234217 -0.854458  0.000000
#         :              :     :      :         :         :
#   }
# }
#
# You can omit the atom charge in your molecule definition.
# (Partial charges will be assigned later according to the force field rules.)
# Responsibility for choosing the atom types (eg "@atom:88", "@atom:89") falls
# on the user.  You must select the type of each atom in the molecule carefully
# by looking at the description in the "In Charges" section of this file
# (see below), and looking for a reasonable match. If your simulation is
# non-neutral, or moltemplate complains that you have missing bond, angle, or
# dihedral types, this means at least one of your atom types is incorrect.


OPLSAA {

  # Below we will use lammps "set" command to assign atom charges
  # by atom type.  https://docs.lammps.org/set.html

  # NOTE1: the commented blocks that you'll find are copied as found in the 
  #   original FF-file, so they don't respect the format/syntax used here
  #   (I thought some of them could be useful anyway, so I kept them here)

"""

file_equivalences_header = """

  # ---------- EQUIVALENCE CATEGORIES for bonded interaction lookup ----------
  #   Each type of atom has a separate ID used for looking up bond parameters
  #   and a separate ID for looking up 3-body angle interaction parameters
  #   and a separate ID for looking up 4-body dihedral interaction parameters
  #   and a separate ID for looking up 4-body improper interaction parameters
  #   The complete @atom type name includes ALL of these ID numbers.  There's
  #   no need to force the end-user to type the complete name of each atom.
  #   The "replace" command used below informs moltemplate that the short
  #   @atom names we have been using above are equivalent to the complete
  #   @atom names used below:

"""

nb_header = """

  # --------------- Non-Bonded interactions: ---------------------
  # https://docs.lammps.org/pair_lj.html
  # Syntax:
  # pair_coeff    AtomType1    AtomType2   parameters...

"""

bond_header = """

  # ------- Bond Interactions: -------
  # https://docs.lammps.org/bond_harmonic.html
  # Syntax:  
  # bond_coeff BondTypeName  parameters...

"""

angle_header = """

  # ------- Angle Interactions: -------
  # https://docs.lammps.org/angle_harmonic.html
  # Syntax:  
  # angle_coeff AngleTypeName  parameters...

"""

dihedral_header = """

  # ----------- Dihedral Interactions: ------------
  # https://docs.lammps.org/dihedral_opls.html
  # Syntax:
  # dihedral_coeff DihedralTypeName  parameters...

"""

improper_header = """

  # ---------- Improper Interactions: ----------
  # https://docs.lammps.org/dihedral_opls.html
  # https://docs.lammps.org/improper_cvff.html
  # https://docs.lammps.org/improper_harmonic.html
  # NOTE: impropers are a WIP, currently implemented as improper_cvff
  # Syntax:
  # improper_coeff ImproperTypeName  parameters

"""

closing_stuff = """

  # LAMMPS supports many different kinds of bonded and non-bonded
  # interactions which can be selected at run time.  Eventually
  # we must inform LAMMPS which of them we will need.  We specify
  # this in the "In Init" section: 

  write_once("In Init") {
    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style opls
    improper_style cvff
    # NOTE: in the original oplsaa.lt file the pair style was
    #   lj/cut/coul/long 11.0 11.0
    # but with an accompanying note stating that OPLSAA/M (2015) 
    # uses a different pair style, the one used here
    # (as I trusted the original author)
    pair_style lj/charmm/coul/long 9.0 11.0
    pair_modify mix geometric
    special_bonds lj/coul 0.0 0.0 0.5
    kspace_style pppm 0.0001
  } #end of init parameters

"""

data_from_atm_num = {
    0: {"element": "XX", "mass": "0.00000000000000001"},
    1: {"element": "H", "mass": "1.008"},
    2: {"element": "He", "mass": "4.003"},
    3: {"element": "Li", "mass": "6.941"},
    4: {"element": "Be", "mass": "9.012"},
    5: {"element": "B", "mass": "10.811"},
    6: {"element": "C", "mass": "12.011"},
    7: {"element": "N", "mass": "14.007"},
    8: {"element": "O", "mass": "15.999"},
    9: {"element": "F", "mass": "18.998"},
    10: {"element": "Ne", "mass": "20.179"},
    11: {"element": "Na", "mass": "22.990"},
    12: {"element": "Mg", "mass": "24.305"},
    13: {"element": "Al", "mass": "26.982"},
    14: {"element": "Si", "mass": "28.086"},
    15: {"element": "P", "mass": "30.974"},
    16: {"element": "S", "mass": "32.065"},
    17: {"element": "Cl", "mass": "35.453"},
    18: {"element": "Ar", "mass": "39.948"},
    19: {"element": "K", "mass": "39.098"},
    20: {"element": "Ca", "mass": "40.078"},
    30: {"element": "Zn", "mass": "65.377"},
    35: {"element": "Br", "mass": "79.904"},
    36: {"element": "Kr", "mass": "83.798"},
    37: {"element": "Rb", "mass": "85.468"},
    38: {"element": "Sr", "mass": "87.620"},
    53: {"element": "I", "mass": "126.905"},
    54: {"element": "Xe", "mass": "131.293"},
    55: {"element": "Cs", "mass": "132.905"},
    56: {"element": "Ba", "mass": "137.327"},
    57: {"element": "La", "mass": "138.905"},
    60: {"element": "Nd", "mass": "144.242"},
    63: {"element": "Eu", "mass": "151.964"},
    64: {"element": "Gd", "mass": "157.25"},
    70: {"element": "Yb", "mass": "173.04"},
    89: {"element": "Ac", "mass": "227.000"},
    90: {"element": "Th", "mass": "232.038"},
    92: {"element": "U", "mass": "238.029"},
    95: {"element": "Am", "mass": "243.000"},
    99: {"element": "DM", "mass": "1.000"},
}
