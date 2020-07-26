#!/bin/bash

# Author: Otello M Roscioni  (roscioni at gmail)
# License: MIT License  (See LICENSE.md) 
# Copyright (c) 2020

# Usage: molc.sh /path/to/In\ Settings

# temporary directory
temp=output_molc.$$

# Set Internal Field Separator to "newline".
IFS=$'\n'

# 0. Commands which may be present in the Settings and must be defined *before* the pair styles.
grep --color=none -e 'replicate' "$1"

# 1. Print the masses.
grep --color=none "mass " "$1"
echo

# 2.1 Check the boundary keyword.
#qstyle=$(awk '{if(tolower($2)~/hybrid\/overlay/){printf "%s",$8; exit}}' "$1") # OLD way: use the Coulomb style specified in the settings section.
bstyle=$(awk '{if(tolower($0)~/boundary/){printf "%s %s %s",tolower($2),tolower($3),tolower($4); exit}}' "$2")
kstyle_mod=""
if [ "$bstyle" == "f f f" ]
then
  qstyle="coul/cut/offcentre"
elif [ "$bstyle" == "p p f" ]
then
  qstyle="coul/long/offcentre"
  kstyle_mod="kspace_modify slab 3.0"
elif [ "$bstyle" == "p p p" ]
then
  qstyle="coul/long/offcentre"
else
  echo "MOLC SYNTAX ERROR: Incorrect boundary conditions." >&2
  echo "Check the file In Init or your input LT script." >&2
  echo >&2
  echo "OFFENDING LINE: " $bstyle >&2
  exit 0
fi

# 2.2 Merge the pair_styles hybrid/overlay gayberne plus coul/long/offcentre.
style=$(awk '{if(tolower($2)~/hybrid\/overlay/){printf "%s %s %s %g %g %g %g %s %g",$1,tolower($2),tolower($3),$4,$5,$6,$7,tolower($8),$9; exit}}' "$1")
cutoff=$(awk '{if(tolower($2)~/hybrid\/overlay/){printf "%g",$9; exit}}' "$1")
rm -f $temp
# don't break the line at spaces.
for string in $(awk '{if(tolower($2)~/hybrid\/overlay/){print}}' "$1"); do
  
  # Consistency test.
  test=$(echo $string | awk '{printf "%s %s %s %g %g %g %g %s %g",$1,tolower($2),tolower($3),$4,$5,$6,$7,tolower($8),$9}')
  if [ "$test" != "$style" ]; then 
  echo "MOLC SYNTAX ERROR: pair_style differs between atom types:" >&2
  echo $style >&2
  echo $test >&2
  echo "OFFENDING LINE: " $string >&2
  exit 0
  fi
  
  # Store the offcentre charges.
  echo $string |\
  awk '{for (i=11; i<=NF; i++){printf("%s ",$i)}}' >> $temp
done
nq=$(wc -w $temp|awk '{print $1/5}') # total number of offcentre charges.

# 2.3 Change the Coulomb style according to qstyle.
out=$(echo $style | awk -v qs=$qstyle '{printf "%s %s %s %g %g %g %g %s %g",$1,$2,$3,$4,$5,$6,$7,qs,$9}')
awk -v a=$out -v b=$nq '{printf "%s %i %s\n",a,b,$0}' $temp
echo

# 3. Long-range electrostatics (if required).
if [ "$qstyle" == "coul/long/offcentre" ]; then
  kstyle=$(awk '{if(tolower($1)~/kspace_style/){printf "%s %s",$2,$3; exit}}'  "$1")
  if [ -z "$kstyle" ]; then
    echo "MOLC SYNTAX ERROR: kspace_style missing" >&2
    echo "Inconsistent bounday conditions." >&2
    exit 0
  fi
  awk -v a=$kstyle -v b=$nq '{printf "kspace_style %s %i %s\n",a,b,$0}' $temp
  if [ -n "$kstyle_mod" ]; then
    echo $kstyle_mod
  fi
fi
echo

# 4. Print the pair coefficients and their explicit combinations.
echo "# Pair coefficients: Coulomb and Gay-Berne parameters."
if [ "$qstyle" == "coul/long/offcentre" ]; then
  echo "pair_coeff * * $qstyle"
else
  echo "pair_coeff * * $qstyle $cutoff"
fi
# Store the GB coefficients for the same atom type into an array.
gb_coeff=( $(awk '{if(tolower($0)~/pair_coeff/ && tolower($0)~/gayberne/ && $2~$3){print}}' "$1") )
# Print the mixed combintations of GB coefficients.
iend=$(bc<<<${#gb_coeff[*]}-2)
jend=$(bc<<<${#gb_coeff[*]}-1)
for i in $(seq 0 $iend); do
  jstart=$(bc<<<$i+1)
  for j in $(seq $jstart $jend); do
    echo ${gb_coeff[$i]} ${gb_coeff[$j]} |\
    awk '{printf "pair_coeff %s %s gayberne %.6g %.6g %s %s %s %s %s %s %s\n",$2,$15,sqrt($5*$18),($6+$19)/2,$7,$8,$9,$20,$21,$22,$13}'
  done
done
# Print the GB coefficients for the same atom type.
awk '{if(tolower($0)~/pair_coeff/ && tolower($0)~/gayberne/ && $2~$3){print}}' "$1" 
echo

# 5. Bonded potential (if required).
grep --color=none -e 'bond_style' "$1" |head -1
grep --color=none -e 'special_bonds' "$1" |head -1
echo
grep --color=none -e 'bond_coeff' "$1"
echo

# 6. Extra Stuff which may be present in the Settings.
grep --color=none -e 'group\|neigh\' "$1"

# final clean-up.
rm -f $temp
