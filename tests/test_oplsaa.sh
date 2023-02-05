#!/usr/bin/env bash

test_oplsaa() {
  cd tests/
    cp -r ../examples/all_atom/force_field_OPLSAA/ethylene+benzene .
    cd ethylene+benzene/
      bash README_setup.sh
      assertTrue "system.data file not created" "[ -s system.data ]"
      NUM_IMPROPERS=`grep impropers system.data | awk '{print $1}'`
      assertTrue "system.data missing impropers" "[ $NUM_IMPROPERS -gt 0 ]"
      cleanup_moltemplate.sh
      NUM_IMPROPERS=`grep impropers system.data | awk '{print $1}'`
      assertTrue "cleanup_moltemplate.sh failed: system.data missing impropers" "[ $NUM_IMPROPERS -gt 0 ]"
      NUM_IMPROPER_TYPES=`grep "improper types" system.data | awk '{print $1}'`
      assertTrue "cleanup_moltemplate.sh failed: wrong number of improper types after cleanup" "[ $NUM_IMPROPER_TYPES -eq 2 ]"
      assertTrue "cleanup_moltemplate.sh failed: system.in.charges file not created" "[ -s system.in.charges ]"
    cd ../
    rm -rf ethylene+benzene/
  cd ../
}

. tests/shunit2/shunit2
