#!/usr/bin/env bash

test_oplsaa() {
  cd tests/
    cd test_molc_files/
      moltemplate.sh -atomstyle "hybrid molecular ellipsoid" -molc system.lt
      assertTrue "system.data file not created" "[ -s system.data ]"
      assertTrue "system.in.settings file not created" "[ -s system.in.settings ]"
    cd ../
  cd ../
}

. shunit2/shunit2
