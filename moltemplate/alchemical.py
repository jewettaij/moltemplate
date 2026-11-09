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
g_date_str = '2026-11-09'
g_version_str = '0.2.0'
import sys
import re
import os
try:
    from .lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid
    from .ttree_lex import SplitQuotedString
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid
    from ttree_lex import SplitQuotedString


def load_type_mappings(in_alchemical_path):
    """
    The content of the mapping file is based on the LAMMPS set atom command, as in:
    set atom $atom:mols[2]/c1  type @atom:OPLSAA/157
    """
    mappings = []  # New atom type mappings
    shape_mappings = []  # New atom shape mappings
    with open(in_alchemical_path, 'r', encoding='utf-8') as f:
        for line_orig in f:
            line = line_orig.rstrip('\n')
            comment = ''
            if '#' in line_orig:
                ic = line.find('#')
                line = line_orig[:ic]
                comment = ' ' + line_orig[ic:].rstrip('\n')
            tokens = SplitQuotedString(line.strip())
            if len(tokens) == 0:
                continue
            # Process the "set atom" commands
            if tokens[0].lower() == 'set' and tokens[1].lower() == 'atom':
                number_tokens = len(tokens)
                # add (key, val) tuples to mappings:
                if tokens[3].lower() == 'type' and number_tokens >= 5:
                    key = tokens[2]
                    val = tokens[4]
                    mappings.append((key, val))
                elif tokens[3].lower() == 'shape' and number_tokens >= 7:
                    key = tokens[2]
                    val = f"{tokens[4]} {tokens[5]} {tokens[6]}"
                    shape_mappings.append((key, val))

    return mappings, shape_mappings



def make_regex_from_key(key):
    """
    Create a compiled regular expression that can recognize strings like
    "$atom:mols[0]/C1", "$atom:mols[1]/C1", or "$atom:mols[2]/C2", ...
    given a ``key`` argument which can take values like:
                  "$/atom:mols[2]/C1", "$/atom:mols[*]/C1"
                  "$/atom:mols[0][*]/C1", "$/atom:mols[*][*]/C1", etc...
    Replace every "[*]" token with a regex that matches a numeric index
    while escaping all other characters.

    NOTE: Currently as of 2026-11-09, the moltemplate.sh script does not allow
    wildcards (e.g. "[*]") to appear in "set" commands.  (You will get a
    "class/variable "mols[*]" not found" error if you try).  So you cannot put
    "set atom $atom:mols[*]/c1 ..." inside a write("Postprocess") command.

    Consequently, the following simpler version of this function
    would work just as well:
        def make_regex_from_key(atom: str) -> re.Pattern:
            escaped = re.escape(atom.lstrip("$"))
            return re.compile(r"\s*\$" + escaped)
    ...But I will leave the complicated version of this function in place in
    case somebody eventually wants to fix this limitation in moltemplate.sh.
    """
    i = 0
    if key.find("$/atom:") != -1:
        i = len("$/atom:")
    atom_name_pattern = key[i:]   # e.g. "mols[2]/C1", or "mols[*][*]/C1"
    parts = atom_name_pattern.split('[*]')
    if len(parts) == 1:
        key_regex = re.escape(atom_name_pattern)
    else:
        pieces = []
        for i, part in enumerate(parts):
            pieces.append(re.escape(part))
            if i != len(parts) - 1:
                pieces.append(r'\[\d+\]')
        key_regex = ''.join(pieces)
    pattern = r'^\$/atom:' + key_regex + r'\b'  # ?
    return re.compile(pattern)


def update_data_atoms(data_atoms_path, type, mappings):
    # compile regexes
    compiled = [(make_regex_from_key(k), v) for k, v in mappings]

    out_lines = []
    num_changed = 0
    with open(data_atoms_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip('\n')
            for rx, val in compiled:
                if rx.search(line):
                    parts = line.split()
                    if len(parts) >= 3:
                        parts[type] = val
                        line = ' '.join(parts)
                        num_changed += 1
                    break
            out_lines.append(line)

    if num_changed:
        with open(data_atoms_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines) + '\n')

    return num_changed


def update_ellipsoids(ellipsoid_path, shape_mappings):
    # compile regexes
    compiled = [(make_regex_from_key(k), v) for k, v in shape_mappings]

    out_lines = []
    num_changed = 0
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
                        num_changed += 1
                    break
            out_lines.append(line)

    if num_changed:
        with open(ellipsoid_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines) + '\n')

    return num_changed
 

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
    num_changed = update_data_atoms(data_path, atom_type, mappings)
    print(f'Alchemical transformation...\nUpdated {num_changed} lines in {data_path}')

    # Update ellipsoids section, if present
    if os.path.isfile(ellipsoid_path) and shape_mappings:
        num_changed = update_ellipsoids(ellipsoid_path, shape_mappings)
        print(f'Updated {num_changed} lines in {ellipsoid_path}')
    
    print("done\n")


if __name__ == '__main__':
    main()

