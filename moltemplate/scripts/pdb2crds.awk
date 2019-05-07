#!/usr/bin/awk -f

# Extracts the x,y,z coordinates from ALL ATOM and HETATM records.
#
# Usage:
#    pdb2crds.awk < input_file.pdb > output_file.txt

/^ATOM  |^HETATM/{print substr($0,31,8)" "substr($0,39,8)" "substr($0,47,8)}
