# Build a new DATA file which contains the topology from the original 
# system.data file, and the coordinates from the latest simulation result
# (which should be stored in the "system_linked_length=1700nm.data" file).

# You can do this by editing the text files by hand,
# extracting the text in the "Atoms" section from the new data file
# and pasting it into the old data file (while also updating the number of
# bond types).  However the following commands do this automatically
#
# (Note "extract_lammps_data.py" is a program that comes with moltemplate.)

TARGET_FILE="system_length=1700nm.data"

extract_lammps_data.py Header < system.data > old_header.tmp
extract_lammps_data.py Header < system_linked_length=1700nm.data >new_header.tmp
extract_lammps_data.py Masses < system.data > old_masses.tmp
extract_lammps_data.py Atoms < system_linked_length=1700nm.data > new_atoms.tmp
extract_lammps_data.py -n Header Masses Atoms < system.data > old_after_atoms.tmp

cat old_header.tmp  > "$TARGET_FILE"
grep "bond types" new_header.tmp >> "$TARGET_FILE"
grep " xlo " new_header.tmp >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"
echo "Masses   # full" >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"
cat old_masses.tmp >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"
echo "Atoms   # full" >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"
cat new_atoms.tmp >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"
cat old_after_atoms.tmp >> "$TARGET_FILE"
echo "" >> "$TARGET_FILE"

rm -f old*.tmp new*.tmp

