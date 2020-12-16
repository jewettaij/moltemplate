#!/bin/sh

SYNTAX_MSG=$(cat <<EOF
Typical Usage:

amberfrcmod2lt.sh DABPA_AC.frcmod GAFF > dabpa_ac_frcmod.lt

(Then include the line "import dabpa_ac_frcmod.lt" at the beginning of the LT
 describing the molecule(s) that needs this frcmod.  For example

import "gaff.lt"
import "dabpa_ac_frcmod.lt"
DABPA inherits GAFF {
  ...define molecule here...
}  

You may need to manually split the .frcmod file and run these scripts instead:
amberparm_pair_to_lt.py, amberparm_bond_to_lt.py, amberparm_angle_to_lt.py...)
(Be sure that all of these .py files are in your PATH as well.)

EOF
)

if [ "$#" != "2" ]; then
    echo "${SYNTAX_MSG}" >&2
    echo "" >&2
    echo "Error: This script requires two arguments," >&2
    echo "       1) the name of the amber parm file to be converted (eg \"gaff.dat\")" >&2
    echo "       2) the name of the moltemplate object to be created (eg \"GAFF\")" >&2
    echo "          (This may include the \"inherits\" keyword and parent classes.)" >&2
    exit 1
fi

MOLTEMPLATE_USAGE_MSG=$(cat <<EOF
#    Background information and usage explanation:
# This file contanis a list of atom types and rules for generating bonded
# interactions between these atoms (hopefully) according to AMBER conventions.
# By using the atom types shown below in your own molecules, bonds and angular
# interactions will be automatically generated.
# AMBER (GAFF) force-field parameters will also be assigned to each angle
# interaction (according to these atom types).
# One way to apply the GAFF force field to a particular type of molecule, is
# to use the "inherits" keyword when you define that molecule.  For example:
# import("gaff.lt")
# MoleculeType inherits GAFF {
#   write_once("Data Atoms") {
#     \$atom:C1 \$mol:... @atom:cx 0.0 4.183 3.194 13.285
#     \$atom:C2 \$mol:... @atom:cx 0.0 4.291 4.618 13.382
#        :       :         :
#   }
# }
#(See "Inheritance" and "short names vs. full names" in the moltemplate manual.)
EOF
)
# (Note that the full name of the atom type in this example is "@atom:/GAFF/cx"
#  You can always refer to atom types this way as well.  Using "inherits GAFF"
#  allows you to use more conventient "@atom:cx" shorthand notation instead.)

echo "####################################################################"
echo "# To use this, LAMMPS currently must be compiled with the USER-MISC package."
echo "# (Type \"make yes-user-misc\" into the shell before compiling LAMMPS.)"
echo "####################################################################"
echo "#    This moltemplate (LT) file was generated automatically using"
echo "# amberfrcmod2lt.sh $1 $2"
echo "####################################################################"
echo "$MOLTEMPLATE_USAGE_MSG"
echo "####################################################################"
echo "#    Moltemplate can not assign atom charge.  You must assign atomic"
echo "# charges yourself.  (Moltemplate is only a simple text manipulation tool.)"
echo "# You can do this afterwards using commands like \"set atom 70 charge -0.212\""
echo "# For details, see http://lammps.sandia.gov/doc/set.html)"
echo "####################################################################"


if ! which ./amberparm_atomdescr_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_atomdescr_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_mass_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_mass_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_pair_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_pair_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_bond_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_bond_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_angle_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_angle_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_dihedral_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_dihedral_to_lt.py\" not found.\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi
if ! which ./amberparm_improper_to_lt.py > /dev/null; then
    echo "\nError: \"amberparm_improper_to_lt.py\" not found. (Update your PATH?)\n" >&2
    echo "       (Try running this script from the directory containing amberparm2lt.sh)" >&2
    exit 2
fi


#FRCMOD_FILE='gaff.dat'
FRCMOD_FILE=$1

# sections are separated by blank lines
# some sections have comment lines at the beginning

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="MASS")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.mass

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="BOND")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.bond

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="ANGLE")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.angle

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="DIHE")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.dihedral

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="IMPROPER")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.improper

tail -n +2 < "$FRCMOD_FILE" | awk 'BEGIN{ignore=1} {if ((NF==1)&&($1=="NONBON")) {ignore=0} else {if (NF==1) {ignore=1} if (!ignore) {if (NF>0) print $0}}}' > ${FRCMOD_FILE}.pair


./amberparm_atomdescr_to_lt.py < "${PARM_FILE}.mass" > "${PARM_FILE}.atomdescr.lt"
./amberparm_mass_to_lt.py < "${FRCMOD_FILE}.mass" > "${FRCMOD_FILE}.mass.lt"
./amberparm_pair_to_lt.py < "${FRCMOD_FILE}.pair" > "${FRCMOD_FILE}.pair.lt"
./amberparm_bond_to_lt.py < "${FRCMOD_FILE}.bond" > "${FRCMOD_FILE}.bond.lt"
./amberparm_angle_to_lt.py < "${FRCMOD_FILE}.angle" > "${FRCMOD_FILE}.angle.lt"
./amberparm_dihedral_to_lt.py \
     < "${FRCMOD_FILE}.dihedral" > "${FRCMOD_FILE}.dihedral.lt"
./amberparm_improper_to_lt.py \
     < "${FRCMOD_FILE}.improper" > "${FRCMOD_FILE}.improper.lt"


echo "#"
echo "# --- Description of atom types ---"
echo "#"
cat "${PARM_FILE}.atomdescr.lt"
echo ""
echo ""
echo ""


echo "$2 {"
echo ""
echo "  # ----------------------------------------------------------------------"
#echo "  # This file was automatically generated by \"common/amber/amberparm2lt.sh\""
echo "  # The basic atom nomenclature and conventions are explained here:"
echo "  #   http://ambermd.org/antechamber/gaff.pdf"
echo "  # For reference, the original gaff.dat file and format documentation are here:"
echo "  #   http://ambermd.org/AmberTools-get.html"
echo "  #   http://ambermd.org/formats.html#parm.dat"
echo "  # ----------------------------------------------------------------------"
echo ""

cat "$FRCMOD_FILE.mass.lt" \
    "$FRCMOD_FILE.pair.lt" \
    "$FRCMOD_FILE.bond.lt" \
    "$FRCMOD_FILE.angle.lt" \
    "$FRCMOD_FILE.dihedral.lt" \
    "$FRCMOD_FILE.improper.lt"

# Commenting out the following (should be already defined in main FF file.)
#
#AMBER_STYLES_INIT=$(cat <<EOF
#
#  write_once("In Init") {
#    # Default styles and settings for AMBER based force-fields:
#    units           real
#    atom_style      full
#    bond_style      hybrid harmonic
#    angle_style     hybrid harmonic
#    dihedral_style  hybrid fourier
#    improper_style  hybrid cvff
#    pair_style      hybrid lj/charmm/coul/long 9.0 10.0 10.0
#    kspace_style    pppm 0.0001
#
#    # NOTE: If you do not want to use long-range coulombic forces,
#    #       comment out the two lines above and uncomment this line:
#    # pair_style      hybrid lj/charmm/coul/charmm 9.0 10.0
#
#    pair_modify     mix arithmetic
#    special_bonds   amber
#  }
#EOF
#)
#
#echo "$AMBER_STYLES_INIT"

echo ""
echo "}"
echo ""
echo ""

# Delete the temporary files we created earlier
rm -f "$FRCMOD_FILE.mass.lt"
rm -f "$FRCMOD_FILE.pair.lt"
rm -f "$FRCMOD_FILE.bond.lt"
rm -f "$FRCMOD_FILE.angle.lt"
rm -f "$FRCMOD_FILE.dihedral.lt"
rm -f "$FRCMOD_FILE.improper.lt"
