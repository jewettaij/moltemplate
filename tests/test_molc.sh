#!/usr/bin/env bash

test_molc() {
  cd tests/
    cd test_molc_files/
    
      # Legagy test.
      moltemplate.sh -atomstyle "hybrid molecular ellipsoid" -molc system.lt
      assertTrue "system.data file not created" "[ -s system.data ]"
      assertTrue "system.in.settings file not created" "[ -s system.in.settings ]"
      
      # Latest version (pair_style molc/cut and pair_style molc/long)
      moltemplate.sh -atomstyle "hybrid molecular ellipsoid" -overlay-bonds -dump water_start.dump sample.lt
      assertTrue "sample.data file not created" "[ -s sample.data ]"
      assertTrue "sample.in.settings file not created" "[ -s sample.in.settings ]"
    cd ../
  cd ../
}

. tests/shunit2/shunit2
