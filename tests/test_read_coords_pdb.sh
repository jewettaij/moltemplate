#!/usr/bin/env bash

test_read_coords_pdb() {
  cd tests/
    cp -r ../examples/file_conversion_examples/read_PDB_file_examples/waterSPCE_from_PDBfile .
    cd waterSPCE_from_PDBfile/
      bash README_setup.sh
      assertTrue "system.data file not created" "[ -s system.data ]"
      NUM_ATOMS=`grep atoms system.data | awk '{print $1}'`
      assertTrue "system.data missing atoms" "[ $NUM_ATOMS -gt 0 ]"
      NUM_BONDS=`grep bonds system.data | awk '{print $1}'`
      assertTrue "system.data missing bonds" "[ $NUM_BONDS -gt 0 ]"
      NUM_ANGLES=`grep angles system.data | awk '{print $1}'`
      assertTrue "system.data missing angles" "[ $NUM_ANGLES -gt 0 ]"
      
      # Extract the coordinates of the first atom of the DATA file (system.data)
      # and confirm that they actually match the coordinates from the PDB file.
      # (We only do this for one of the lines.  I chose line #4 arbitrarily.)
      CRDS_DATA_LINE4=`extract_lammps_data.py Atoms < system.data | sort -g -k 1 | awk '{print $5" "$6" "$7}' | awk '{if (NR==4) {print $0}}'`
      CRDS_PDB_LINE4=`awk '/^ATOM  |^HETATM/{print substr($0,31,8)" "substr($0,39,8)" "substr($0,47,8)}' < moltemplate_files/solvate.pdb  | awk '{if (NR==4) {print $0}}'`
      Xdata=`echo $CRDS_DATA_LINE4 | awk '{print $1}'`
      Ydata=`echo $CRDS_DATA_LINE4 | awk '{print $2}'`
      Zdata=`echo $CRDS_DATA_LINE4 | awk '{print $3}'`
      Xpdb=`echo $CRDS_PDB_LINE4 | awk '{print $1}'`
      Ypdb=`echo $CRDS_PDB_LINE4 | awk '{print $2}'`
      Zpdb=`echo $CRDS_PDB_LINE4 | awk '{print $3}'`
      assertTrue "PDB coordinates do not match coordinates from system.data file created by moltemplate. (See line 4)" "[ $Xdata = $Xpdb ] && [ $Ydata = $Ypdb ] && [ $Zdata = $Zpdb ]"
    cd ../
    rm -rf waterSPCE_from_PDBfile/
  cd ../
}

. shunit2/shunit2
