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
    cd ../
    rm -rf waterSPCE_from_PDBfile/
  cd ../
}

. shunit2/shunit2
