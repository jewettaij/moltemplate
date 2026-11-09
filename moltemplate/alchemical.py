#!/usr/bin/env python3

# Author: Otello M Roscioni (https://orcid.org/0000-0001-7815-6636)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2025
"""
Perform an alchemical transformation.

Before the remove_duplicate_atoms.py script is invoked in Moltemplate.sh,
look for specific atom
"""
g_program_name = __file__.split('/')[-1]  # ='alchemical.py'
g_date_str = '2026-01-16'
g_version_str = '0.1.1'

import sys
import re
import os
try:
    from .lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid
except Exception:
    from lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid

def load_type_mappings(in_alchemical_path):
    """
    The content of the mapping file is based on the LAMMPS set atom command, as in:
    
    set atom mol1[2]/C1 type @atom:OPLSAA/99
    """
    mappings = []  # New atom type mappings
    shape_mappings = []  # New atom shape mappings
    with open(in_alchemical_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            # Process the "set atom" commands
            if parts[0].lower() == 'set' and parts[1].lower() == 'atom':
                number_parts = len(parts)
                # add (key, val) tuples to mappings:
                if parts[3].lower() == 'type' and number_parts >= 5:
                    key = parts[2]
                    val = parts[4]
                    mappings.append((key, val))
                elif parts[3].lower() == 'shape' and number_parts >= 7:
                    key = parts[2]
                    val = f"{parts[4]} {parts[5]} {parts[6]}"
                    shape_mappings.append((key, val))

    return mappings, shape_mappings


def make_regex_from_key(key):
    # Examples of keys: mol1[2]/C1, mol1[*]/C1, mol1[0][*]/C1, mol1[*][*]/C1, etc.
    # Replace every literal "[*]" token with a regex that matches a numeric index
    # while escaping all other characters.
    parts = key.split('[*]')
    if len(parts) == 1:
        key_regex = re.escape(key)
    else:
        pieces = []
        for i, part in enumerate(parts):
            pieces.append(re.escape(part))
            if i != len(parts) - 1:
                pieces.append(r'\[\d+\]')
        key_regex = ''.join(pieces)

    pattern = r'^' + key_regex + r'\b'
    return re.compile(pattern)


def update_data_atoms(data_atoms_path, type, mappings):
    # compile regexes
    compiled = [(make_regex_from_key(k), v) for k, v in mappings]

    out_lines = []
    changed = 0
    with open(data_atoms_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip('\n')
            for rx, val in compiled:
                if rx.search(line):
                    parts = line.split()
                    if len(parts) >= 3:
                        parts[type] = val
                        line = ' '.join(parts)
                        changed += 1
                    break
            out_lines.append(line)

    if changed:
        with open(data_atoms_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines) + '\n')

    return changed


def update_ellipsoids(ellipsoid_path, shape_mappings):
    # compile regexes
    compiled = [(make_regex_from_key(k), v) for k, v in shape_mappings]

    out_lines = []
    changed = 0
    with open(ellipsoid_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip('\n')
            for rx, val in compiled:
                if rx.search(line):
                    parts = line.split()
                    if len(parts) == 8:
                        parts[1] = val
                        del parts[2:4]
                        line = ' '.join(parts)
                        changed += 1
                    break
            out_lines.append(line)

    if changed:
        with open(ellipsoid_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines) + '\n')

    return changed
 

def main():

    postprocess_path = sys.argv[1]
    data_path = sys.argv[2]
    ellipsoid_path = sys.argv[3]
    ttree_args = sys.argv[4]

    # Get i_atomtype from atomstyle in ttree_args.
    args = str(ttree_args).split()
    atom_style_string = args[args.index("-atomstyle") + 1] if "-atomstyle" in args else "full"

    # Remove surrounding quotes if present
    if (len(atom_style_string) > 2 and
            ((atom_style_string[0] == atom_style_string[-1]) and
            atom_style_string[0] in ('"', "'"))):
        atom_style_string = atom_style_string[1:-1]

    col_names = AtomStyle2ColNames(atom_style_string)
    _, atom_type, _ = ColNames2AidAtypeMolid(col_names)

    # get mappings
    mappings, shape_mappings = load_type_mappings(postprocess_path)
    if not mappings:
        print('Alchemical transformation failed\nNo mappings found in', postprocess_path)
        sys.exit(0)

    # Update data atoms section
    changed = update_data_atoms(data_path, atom_type, mappings)
    print(f'Alchemical transformation...\nUpdated {changed} lines in {data_path}')

    # Update ellipsoids section, if present
    if os.path.isfile(ellipsoid_path) and shape_mappings:
        changed = update_ellipsoids(ellipsoid_path, shape_mappings)
        print(f'Updated {changed} lines in {ellipsoid_path}')
    
    print("done\n")


if __name__ == '__main__':
    main()

