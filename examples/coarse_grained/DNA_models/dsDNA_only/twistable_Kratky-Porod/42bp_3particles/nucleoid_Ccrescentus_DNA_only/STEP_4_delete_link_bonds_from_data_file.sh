# Build a new DATA file which contains the topology from the original 
# system.data file, and the coordinates from the latest simulation result
# (which should be stored in the "system_linked_length=1900nm.data" file).

# You can do this by editing the text files by hand,
# extracting the text in the "Atoms" section from the new data file
# and pasting it into the old data file (while also updating the number of
# bond types).  However the following commands do this automatically
#
# (Note "extract_lammps_data.py" is a program that comes with moltemplate.)

FILE_NEW_TOPO="system.data"
FILE_COORDS="system_linked_length=1900nm.data"
FILE_TARGET="system_length=1900nm.data"

extract_lammps_data.py Header < "$FILE_NEW_TOPO" > old_header.tmp
extract_lammps_data.py Header < "$FILE_COORDS" > new_header.tmp
extract_lammps_data.py Masses < "$FILE_NEW_TOPO" > old_masses.tmp
extract_lammps_data.py Atoms < "$FILE_COORDS" > new_atoms.tmp
extract_lammps_data.py -n Header Masses Atoms < "$FILE_NEW_TOPO" > old_after_atoms.tmp

cat old_header.tmp  > "$FILE_TARGET"
grep "bond types" new_header.tmp >> "$FILE_TARGET"
grep " xlo " new_header.tmp >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"
echo "Masses   # full" >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"
cat old_masses.tmp >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"
echo "Atoms   # full" >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"
cat new_atoms.tmp >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"
cat old_after_atoms.tmp >> "$FILE_TARGET"
echo "" >> "$FILE_TARGET"

rm -f old*.tmp new*.tmp

