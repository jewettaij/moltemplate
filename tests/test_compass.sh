#!/usr/bin/env bash

test_compass() {
  cd tests/
    cp -r ../examples/all_atom/force_field_COMPASS/alkane_chain_single/ .
    cd alkane_chain_single/
      bash README_setup.sh
      assertTrue "system.data file not created" "[ -s system.data ]"
      NUM_DIHEDRALS=`grep dihedrals system.data | awk '{print $1}'`
      assertTrue "system.data missing dihedrals" "[ $NUM_DIHEDRALS -gt 0 ]"
      cleanup_moltemplate.sh
      NUM_DIHEDRALS=`grep dihedrals system.data | awk '{print $1}'`
      assertTrue "cleanup_moltemplate.sh failed: system.data missing dihedrals" "[ $NUM_DIHEDRALS -gt 0 ]"
      assertTrue "cleanup_moltemplate.sh failed: system.in.charges file not created" "[ -s system.in.charges ]"
    cd ../
    rm -rf alkane_chain_single/
  cd ../
}

. shunit2/shunit2
