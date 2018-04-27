#!/usr/bin/env bash

test_ltemplify() {
  cd tests/
  cp -r ../examples/file_conversion_examples/convert_LAMMPS_to_LT_examples/cnad-cnt .
    cd cnad-cnt
      bash README_step1_run_ltemplify.sh
      assertTrue "cnt.lt file not created" "[ -s cnt.lt ]"
      assertTrue "cnad.lt file not created" "[ -s cnad.lt ]"
      NUM_DIH_CNT=`awk "/dihedral:id/{sum+=1} END{print sum}" < cnt.lt`
      assertTrue "cnt.lt missing dihedrals" "[ $NUM_DIH_CNT -gt 0 ]"
      NUM_DIH_CNAD=`awk "/dihedral:id/{sum+=1} END{print sum}" < cnad.lt`
      assertTrue "cnad.lt missing dihedrals" "[ $NUM_DIH_CNAD -gt 0 ]"
      NUM_DIHT_CNT=`awk "/dihedral_coeff/{sum+=1} END{print sum}" < cnt.lt`
      assertTrue "cnt.lt missing dihedral_coeffs" "[ $NUM_DIHT_CNT -gt 0 ]"
      NUM_DIHT_CNAD=`awk "/dihedral_coeff/{sum+=1} END{print sum}" < cnt.lt`
      assertTrue "cnad.lt missing dihedral_coeffs" "[ $NUM_DIHT_CNAD -gt 0 ]"
      bash README_step2_run_moltemplate.sh
      assertTrue "system.data file not created" "[ -s system.data ]"
      NUM_ATOMS=`grep atoms system.data | awk '{print $1}'`
      assertTrue "system.data missing atoms" "[ $NUM_ATOMS -gt 0 ]"
      NUM_BONDS=`grep bonds system.data | awk '{print $1}'`
      assertTrue "system.data missing bonds" "[ $NUM_BONDS -gt 0 ]"
      NUM_ANGLES=`grep angles system.data | awk '{print $1}'`
      assertTrue "system.data missing angles" "[ $NUM_ANGLES -gt 0 ]"
      NUM_DIHEDRALS=`grep dihedrals system.data | awk '{print $1}'`
      assertTrue "system.data missing dihedrals" "[ $NUM_DIHEDRALS -gt 0 ]"
    cd ../
    rm -rf cnad-cnt
  cd ../
}

. shunit2/shunit2
