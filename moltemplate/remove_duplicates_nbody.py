#!/usr/bin/env python3

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2013

"""
   Get rid of lines containing duplicate bonded nbody interactions in the
   corresponding section of a LAMMPS data file (such as bonds, angles,
   dihedrals and impropers).  Duplicate lines which occur later are
   preserved and the earlier lines are erased.
   (This program reads from sys.stdin.  This program does not parse the entire
    data file.  The text from the relevant section of the LAMMPS file should be
   extracted in advance before it is sent to this program.)

"""

import os
import sys
from collections import defaultdict

try:
    from .ttree_lex import SplitQuotedString
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree_lex import SplitQuotedString


def main():
    in_stream = sys.stdin

    interaction_style = ""
    log_warning_filename = ""
    filter_id_str = ""
    filter_type_str = ""
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    if len(sys.argv) > 3:
        interaction_style = sys.argv[2]  # eg. "Dihedral"
        log_warning_filename = sys.argv[3]  # eg. "warning_duplicate_dihedrals.txt"
    if len(sys.argv) > 4:
        filter_id_str = sys.argv[4]  # eg. "bytype"
    if len(sys.argv) > 5:
        filter_type_str = sys.argv[5]  # eg. "__"
    if (len(sys.argv) not in (2,4,5,6)) or (n < 1):
        sys.stderr.write(
            'Error (remove_duplicates_nbody.py): expected 1,3,4, or 5 arguments.\n')
        sys.exit(-1)

    # Keep track of what interaction(s) were assigned to each set of atoms
    atomids2interactions = defaultdict(list)  # Dict[Tuple[str], List[str])
    atomids2interaction = {}  # of type Dict[Tuple[str], str]

    lines = in_stream.readlines()

    # Start at the end of the file and read backwards.
    # If duplicate lines exist, eliminate the ones that occur earlier in the file.
    i = len(lines)
    while i > 0:
        i -= 1
        line_orig = lines[i]
        line = line_orig.rstrip('\n')
        if '#' in line_orig:
            ic = line.find('#')
            line = line_orig[:ic]

        # Split the line into words (tokens) using whitespace delimeters
        tokens = SplitQuotedString(line,
                                   quotes='{',
                                   endquote='}')

        if len(tokens) == 0:
            del lines[i]  # skip blank lines
        elif len(tokens) == 2 + n:
            atom_ids = tuple(tokens[2:2 + n])
            if atom_ids in atomids2interactions:
                # If an interaction already exists between these atoms
                # then delete this one (this line)
                del lines[i]
            else:
                # Then this is the remaining interaction assigned to these atoms
                atomids2interaction[atom_ids] = line
            # Also keep track of all the interactions we did not keep.
            atomids2interactions[atom_ids].append(line)

    for line in lines:
        sys.stdout.write(line)

    # This portion of the code generates warning messages when duplicate
    # interactions were deleted.  But the code is a little confusing because
    # in order to be reported the interactions might need to satisfy a name
    # filter.  Specifically the id_names ($) or type_names (@) might need
    # to contain some string we are looking for.
    if log_warning_filename != "":
        with open(log_warning_filename, 'w') as log_warning:
            # Did the user ask us to report duplicates?
            for atom_ids, interactions in atomids2interactions.items():
                if len(interactions) == 1:
                    continue
                chosen_interaction = atomids2interaction[atom_ids]
                chosen_interaction_tokens = chosen_interaction.split()
                interaction_id = chosen_interaction_tokens[0]  # eg. "$/dihedral:bytype6910"
                interaction_type = chosen_interaction_tokens[1]  # eg. "@/dihedral:OPLSAA/HC_CM_CM_HC"
                # Look for id_filter_string somewhere inside the interaction_id string.
                if filter_id_str != "" and interaction_id.find(filter_id_str) == -1:
                    # If that string is not found, then skip over this interaction.  This is used
                    # to avoid complaining about interactions which lack a certain string.  (For
                    # example, we only complain about interactions which were generated automatically.
                    # Automatically generated interactions have id strings beginning with "bytype".)
                    continue
                # Look for type_filter_string somewhere inside the interaction_type string.
                if filter_type_str != "" and interaction_type.find(filter_type_str) == -1:
                    # If that string is not found, then skip over this interaction.  This is used
                    # to avoid complaining about interactions which lack a certain string.  (For
                    # example, in the OPLSAA force field, we only complain about "__" interactions.)
                    continue
                log_warning.write(
                    f"Duplicate {interaction_style.lower()}s detected involving atoms:\n"
                    "    "
                    + "\n    ".join(atom_ids) + "\n"
                    "  ...of types:\n"
                    f"    "
                    + "\n    ".join([int_str.split()[1] for int_str in interactions]) + "\n"
                    f"  ...but only this {interaction_style.lower()} type was kept:\n"
                    f"    {atomids2interaction[atom_ids].split()[1]}\n"
                    f"  Was this the correct {interaction_style.lower()} type for these atoms?\n"
                    f"  If not, create an explicit {interaction_style.lower()} interaction between those atoms\n"
                    f"  in the \"Data {interaction_style}\" section of your molecule to override this choice.\n"
                    "\n"
                )
            # If the warnings file is empty, there were no warnings.
            # In that case, delete the file:
            if os.stat(log_warning_filename).st_size == 0:
                os.remove(log_warning_filename)
    return


if __name__ == '__main__':
    main()
