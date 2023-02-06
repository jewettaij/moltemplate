#!/usr/bin/env bash

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2013

# cleanup_moltemplate.sh
# This script attempts to remove irrelevant information from LAMMPS 
# input scripts and data files (such as extra, unneeded atom types 
# and force field parameters).
#
# Unfortunately, the system.data and system.in.settings file which are
# created by moltemplate.sh often contain a lot of irrelevant information,
# such as definition of parameters for atom types defined in some force field
# file that the user is using, but not present in the system they are building
#
# In my experience, this extra information appears to be mostly harmless.
# (Loading this information does not seem to slow down LAMMPS significantly.)
# 
# However it can make visualization difficult in VMD.  (Technically, this
# extra information can take up megabytes of space in the system.data
# and system.in.settings files.  Additionally, when you run LAMMPS, an O(n^2)
# sized table is allocated to store the parameters for interactions between 
# every possible pair of atom types (n atom types), and this occupies
# significantly more memory if n is large.  For example, the "oplsaa.lt" file
# and "oplsaa.prm" (TINKER-format) file both define almost 1000 atom types.)
#
# Usage:
#        cleanup_moltemplate [-base BASE_FILE_NAME]
#
#        Invoke this script from a directory containing these files:
#          system.data, system.in.init, system.in.settings, system.in.charges
#        (If your files don't begin with the name "system", you can pass 
#         the "-base" argument to select files with a different base name.)
#
#        Cleanup_moltemplate.sh will modify these files to remove unnecessary
#        atoms and parameters.  (If your files have other names, you must rename
#        them to match moltemplate file name conventions.)
#
# DO NOT USE THIS SCRIPT ON SIMULATIONS CONTAINING MANY-BODY PAIR STYLES,
# DREIDING-STYLE HYDROGEN BONDS, OR SIMS NEEDING NON-STANDARD AUXILIARY FILES.
# (This script relies on ltemplify.py and inherits its limitations.)



if which python3 > /dev/null; then
    PYTHON_COMMAND='python3'
elif which python > /dev/null; then
    PYTHON_COMMAND='python'
elif which python2 > /dev/null; then
    PYTHON_COMMAND='python2'
fi

# Determine the directory in which the python scripts are located.
# (such as ltemplify.py).  It could either be the directory where the script
# file is located, OR it could be the parent of this directory.
SH_SCR_DIR=`dirname "$0"`
PY_SCR_DIR=$SH_SCR_DIR
if [ ! -s "${PY_SCR_DIR}/ltemplify.py" ]; then
    PY_SCR_DIR="$PY_SCR_DIR/.."
fi


BASE_NAME="system"
LTEMPLIFY_ARGS=""
# Store the list ofo arguments in ARGV, and count them in ARGC
ARGC=0
for A in "$@"; do
    A_FIRSTCHAR="$(echo $A| cut -c 1)"
    # (Note to self: this next line only works in bash, not classic sh.)
    if [ "$A_FIRSTCHAR" = "\$" ]; then
        A="\\$A" # put an extra slash in front to prevent expansion later
    fi
    ARGC=$((ARGC+1))
    eval ARGV${ARGC}=\"$A\"
done


i=0
while [ "$i" -lt "$ARGC" ]; do
    i=$((i+1))
    eval A=\${ARGV${i}}

    if [ "$A" = "-base" ]; then
        if [ "$i" -eq "$ARGC" ]; then
            echo "ERROR cleanup_moltemplate.sh: Base file name expected following -base argument" >&2
            exit 7
        fi
        i=$((i+1))
        eval A=\${ARGV${i}}
        FILE_NAME=$A
        # strip off the ".data" suffix (if present)
        BASE_NAME=`basename "$FILE_NAME" ".data"`
        # strip off the ".lmpdat" suffix (if present)
        BASE_NAME=`basename "$BASE_NAME" ".lmpdat"`

    #else:  If the arguments are not understood in this script, then
    #       pass them on to "ltemplify.py"
    else
        A_FIRSTCHAR="$(echo $A| cut -c 1)"

        if [ "$A_FIRSTCHAR" = "\$" ]; then
            A="\\$A" # put an extra slash in front to prevent expansion later
        fi

        if [ -z "$LTEMPLIFY_ARGS" ]; then
            LTEMPLIFY_ARGS="$A"
        else
            LTEMPLIFY_ARGS="${LTEMPLIFY_ARGS} $A"
        fi
        # Check to see if this string ($A) ends in .lt or .LT
        # If so, then set the base name of the output files
        # to equal the base name of the .LT file being read.
        # (Being careful here.
        #  Sometimes the last argument is not the .lt or .LT file.
        #  Sometimes that file appears earlier in the argument list.
        #  I want to supply a default value.)
        #
        #   Note, in bash you can use:
        # if [ "${LAST_ARG/%.lt/}" -neq "$LAST_ARG" ]; then
        #     OUT_FILE_BASE="${LAST_ARG/%.lt/}"
        # But in the original bourn shell (sh), this does not work.
        # Instead we use a hack involving basename and dirname:

        if [ "$A_FIRSTCHAR" != "-" ]; then
            DN=`dirname "$A"`
            if [ "$DN" = "." ]; then
                DN=""
            else
                DN="${DN}/"
            fi

            BN=`basename "$A" .lt`
            if [ "${DN}${BN}" != "$A" ]; then
                OUT_FILE_BASE="$BN"
            else
                BN=`basename "$A" .LT`
                if [ "${DN}${BN}" != "$A" ]; then
                    OUT_FILE_BASE="$BN"
                fi
            fi
        fi
    fi
done



if [ ! -f "${BASE_NAME}.data" ] || [ ! -f "${BASE_NAME}.in.init" ] || [ ! -f "${BASE_NAME}.in.settings" ];
then
    echo "============================ ERROR ==============================" >&2
    echo "The following files must be exist for this script to work" >&2
    echo "  ${BASE_NAME}.data,  ${BASE_NAME}.in.init,  ${BASE_NAME}.in.settings" >&2
    echo "This script assumes that the files you created with moltemplate begin" >&2
    echo "with \"${BASE_NAME}\".  If those files are absent, it means that you did not run" >&2
    echo "  moltemplate.sh ${BASE_NAME}.lt [...]" >&2
    echo "This usually happens if you named your \"${BASE_NAME}.lt\" file to something else." >&2
    echo "" >&2
    if [ "${BASE_NAME}" == "system" ]; then
        echo "If these files begin with a different string (eg \"BASE_NAME\"), then try" >&2
        echo "running cleanup_moltemplate.sh with the \"-base BASE_NAME\" argument." >&2
        echo "(Note that the beginning of all of these file names must agree.)" >&2
    else
        echo "Try changing your main .lt file to \"${BASE_NAME}.lt\"," >&2
        echo "and try running moltemplate.sh again." >&2
    fi
    echo "" >&2
    echo "Alternatively, you can try renaming ALL of the generated files" >&2
    echo "(including the files ending in:" >&2
    echo "\".data\", \".in.init\", \".in.settings\", and \".in.charges\" if present)" >&2
    echo "  to:" >&2
    echo "${BASE_NAME}.data, ${BASE_NAME}.in.init, ${BASE_NAME}.in.settings (${BASE_NAME}.in.charges if present)" >&2
    echo "================================================================" >&2
    exit 1
fi



PATH_TO_DATA_FILE="."

pushd "$PATH_TO_DATA_FILE"

rm -rf new_lt_file_TMP
mkdir new_lt_file_TMP
cd new_lt_file_TMP

  # now run ltemplify.py
  #
  # If any arguments were passed to this script ("$@") insert them here.
  # (An example of a useful argument would be "-ignore-comments".
  #  It will prevent ltemplify.py from trying to infer atom type names
  #  from comments in the "Masses" section of the DATA file.  Such comments
  #  might contain spaces or special characters which we want to avoid.)

  if ! $PYTHON_COMMAND "${PY_SCR_DIR}/ltemplify.py" $LTEMPLIFY_ARGS "../${BASE_NAME}.in".* "../${BASE_NAME}.data" > "${BASE_NAME}.lt"; then
      echo "ERROR: ltemplify.py failed to parse your files." >&2
      exit 1
  fi

  # This creates a new .LT file named ${BASE_NAME}.lt in the local directory.

  # The ltemplify.py script also does not copy the boundary dimensions.
  # We must do this manually.
  # Extract the header of the data file, reverse the order, and read lines
  # until you have 
  # If you did NOT throw away the "Data Boundary" file usually located in
  # "moltemplate_files/output_ttree/Data Boundary"
  # then you can copy that information from this file into ${BASE_NAME}.lt


  # oops. looks like we don't need this after all
  #function _reverse_lines {
  #    # The following function reverses the order of lines in a file:
  #    # (Neither "tac", nor "tail -r" have cross-platform support.)
  #    awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }'
  #}

  echo "" >> "${BASE_NAME}.lt"
  echo "write_once(\"Data Boundary\") {" >> "${BASE_NAME}.lt"
  # Extract the periodic boundary box dimensions from the 
  # end of the header section of the LAMMPS data file:
  $PYTHON_COMMAND "${PY_SCR_DIR}/extract_lammps_data.py" Header < "../${BASE_NAME}.data" | awk '{if (($3=="xlo") && ($4=="xhi")) {xl=$0} if (($3=="ylo") && ($4=="yhi")) {yl=$0} if (($3=="zlo") && ($4=="zhi")) {zl=$0} if (($4=="xy") && ($5=="xz") && ($6=="yz")) {xtr=$0}} END{print xl; print yl; print zl; if (xtr!="") {print xtr}}' >> "${BASE_NAME}.lt"
  echo "}" >> "${BASE_NAME}.lt"
  echo "" >> "${BASE_NAME}.lt"

  # Now, run moltemplate on this new .LT file.
  # Interpret the "${BASE_NAME}.lt literally. Don't check for duplicates("-overlay...")
  if ! "${SH_SCR_DIR}/moltemplate.sh" "${BASE_NAME}.lt" \
                      -overlay-bonds -overlay-angles \
                      -overlay-dihedrals -overlay-impropers; then
      echo "ERROR: cleanup_moltemplate.sh: unable to convert the simplified LT file to LAMMPS files" >&2
      exit 1
  fi
  # This will create: ${BASE_NAME}.data, ${BASE_NAME}.in.init, ${BASE_NAME}.in.settings

  # That's it.  The new ${BASE_NAME}.data and ${BASE_NAME}.in.settings files should
  # be ready to run in LAMMPS.

  # Special case: "set" commands
  # Typically "set type" or "set atom" commands are used to assign atom charge
  # If there is a ${BASE_NAME}.in.charges file, then it contains these commands
  # however the atom type numbers will be wrong, so we must rewrite it.
  # Replace it with the corresponding commands from the ${BASE_NAME}.in.settings
  # (whose atom type numbers are correct)
  if [ -f "../${BASE_NAME}.in.charges" ]; then
    awk '{ if ((NF >= 5) && ($1 == "set") && ($4 == "charge")) print $0}' \
        < "${BASE_NAME}.in.settings" > "${BASE_NAME}.in.charges"
    # There is no need to remove these lines from ""${BASE_NAME}.in.settings",
    # (because there's no harm to invoke the "set" command twice)
    # ...but if you want to do that, try using a command similar to:
    #sed '/set type/,+1 d' < "${BASE_NAME}.in.settings" > "${BASE_NAME}.in.settings.tmp"
    #mv -f "${BASE_NAME}.in.settings.tmp" "${BASE_NAME}.in.settings"
  fi


  # Now move the ${BASE_NAME}.data and ${BASE_NAME}.in.* files to their original location:
  mv -f "${BASE_NAME}.data" "${BASE_NAME}.in".* ../
  cd ../

# Finally, delete all of the temporary files we generated
rm -rf new_lt_file_TMP
popd

