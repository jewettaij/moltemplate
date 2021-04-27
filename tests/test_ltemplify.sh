#!/usr/bin/env bash

test_ltemplify() {
  cd tests/
  cp -r ../examples/file_conversion_examples/convert_LAMMPS_to_LT_examples/cnad-cnt .

    # test for basic ltemplify.py functionality
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

    # test for the ability to infer type names from comments
    cp -r test_ltemplify_files deleteme
    cd deleteme
      ltemplify.py input_script_w_coeffs.in input_data_no_coeffs.data > out.lt
      NUM_PAIR_COEFFS=`awk "/pair_coeff/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks pair_coeffs" "[ $NUM_PAIR_COEFFS -gt 0 ]"
      FOUND_GROUP_COMMAND=`awk "/group gEthylenes/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the group command" "[ $FOUND_GROUP_COMMAND -gt 0 ]"
      FOUND_SHAKE_COMMAND=`awk "/fix fHCbonds gEthylenes shake/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the shake command" "[ $FOUND_SHAKE_COMMAND -gt 0 ]"
      FOUND_SET_COMMAND=`awk "/set type/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the set command" "[ $FOUND_SET_COMMAND -gt 0 ]"
      N_C_ATOM_NAMES=`awk "/atom:C/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output has not enough C atoms" "[ $N_C_ATOM_NAMES -eq 28 ]"
      N_spaces_ATOM_NAMES=`awk "/atom:H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like atom type names with spaces" "[ $N_spaces_ATOM_NAMES -eq 29 ]"
      N_spaces_BOND_NAMES=`awk "/bond:C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like bond type names with spaces" "[ $N_spaces_BOND_NAMES -eq 6 ]"
      N_spaces_ANGLE_NAMES=`awk "/angle:H_C_C spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like angle type names with spaces" "[ $N_spaces_ANGLE_NAMES -eq 5 ]"
      N_spaces_DIHEDRAL_NAMES=`awk "/dihedral:H_C_C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like dihedral type names with spaces" "[ $N_spaces_DIHEDRAL_NAMES -eq 5 ]"
      N_spaces_IMPROPER_NAMES=`awk "/improper:C_H_H_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like improper type names with spaces" "[ $N_spaces_IMPROPER_NAMES -eq 7 ]"
    cd ../
    rm -rf deleteme/

    cp -r test_ltemplify_files deleteme
    cd deleteme
      ltemplify.py input_script_no_coeffs.in input_data_w_coeffs.data > out.lt
      NUM_PAIR_COEFFS=`awk "/pair_coeff/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks pair_coeffs" "[ $NUM_PAIR_COEFFS -gt 0 ]"
      FOUND_GROUP_COMMAND=`awk "/group gEthylenes/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the group command" "[ $FOUND_GROUP_COMMAND -gt 0 ]"
      FOUND_RATTLE_COMMAND=`awk "/fix fHCbonds gEthylenes rattle/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the rattle command" "[ $FOUND_RATTLE_COMMAND -gt 0 ]"
      FOUND_SET_COMMAND=`awk "/set type/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the set command" "[ $FOUND_SET_COMMAND -gt 0 ]"
      N_C_ATOM_NAMES=`awk "/atom:C/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output has not enough C atoms" "[ $N_C_ATOM_NAMES -eq 28 ]"
      N_spaces_ATOM_NAMES=`awk "/atom:H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like atom type names with spaces" "[ $N_spaces_ATOM_NAMES -eq 29 ]"
      N_spaces_BOND_NAMES=`awk "/bond:C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like bond type names with spaces" "[ $N_spaces_BOND_NAMES -eq 6 ]"
      N_spaces_ANGLE_NAMES=`awk "/angle:H_C_C spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like angle type names with spaces" "[ $N_spaces_ANGLE_NAMES -eq 5 ]"
      N_spaces_DIHEDRAL_NAMES=`awk "/dihedral:H_C_C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like dihedral type names with spaces" "[ $N_spaces_DIHEDRAL_NAMES -eq 5 ]"
      N_spaces_IMPROPER_NAMES=`awk "/improper:C_H_H_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like improper type names with spaces" "[ $N_spaces_IMPROPER_NAMES -eq 7 ]"
    cd ../
    rm -rf deleteme/

    cp -r test_ltemplify_files deleteme
    cd deleteme
      ltemplify.py -datacoeffs input_script_no_coeffs.in input_data_w_coeffs.data > out.lt
      FOUND_DATA_PAIR_COEFFS=`awk "/Pair Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data Pair Coeffs section" "[ $FOUND_DATA_PAIR_COEFFS -gt 0 ]"
      FOUND_DATA_PAIRIJ_COEFFS=`awk "/PairIJ Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data PairIJ Coeffs section" "[ $FOUND_DATA_PAIRIJ_COEFFS -gt 0 ]"
      FOUND_DATA_BOND_COEFFS=`awk "/Bond Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data Bond Coeffs section" "[ $FOUND_DATA_BOND_COEFFS -gt 0 ]"
      FOUND_DATA_ANGLE_COEFFS=`awk "/Angle Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data Angle Coeffs section" "[ $FOUND_DATA_ANGLE_COEFFS -gt 0 ]"
      FOUND_DATA_DIHEDRAL_COEFFS=`awk "/Dihedral Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data Dihedral Coeffs section" "[ $FOUND_DATA_DIHEDRAL_COEFFS -gt 0 ]"
      FOUND_DATA_IMPROPER_COEFFS=`awk "/Improper Coeffs/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output lacks a Data Improper Coeffs section" "[ $FOUND_DATA_IMPROPER_COEFFS -gt 0 ]"
      FOUND_GROUP_COMMAND=`awk "/group gEthylenes/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the group command" "[ $FOUND_GROUP_COMMAND -gt 0 ]"
      FOUND_RATTLE_COMMAND=`awk "/fix fHCbonds gEthylenes rattle/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the rattle command" "[ $FOUND_RATTLE_COMMAND -gt 0 ]"
      FOUND_SET_COMMAND=`awk "/set type/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py failed to process the set command" "[ $FOUND_SET_COMMAND -gt 0 ]"
      N_C_ATOM_NAMES=`awk "/atom:C/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py output has not enough C atoms" "[ $N_C_ATOM_NAMES -eq 28 ]"
      N_spaces_ATOM_NAMES=`awk "/atom:H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like atom type names with spaces" "[ $N_spaces_ATOM_NAMES -eq 29 ]"
      N_spaces_BOND_NAMES=`awk "/bond:C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like bond type names with spaces" "[ $N_spaces_BOND_NAMES -eq 6 ]"
      N_spaces_ANGLE_NAMES=`awk "/angle:H_C_C spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like angle type names with spaces" "[ $N_spaces_ANGLE_NAMES -eq 5 ]"
      N_spaces_DIHEDRAL_NAMES=`awk "/dihedral:H_C_C_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like dihedral type names with spaces" "[ $N_spaces_DIHEDRAL_NAMES -eq 5 ]"
      N_spaces_IMPROPER_NAMES=`awk "/improper:C_H_H_H spaces/{sum+=1} END{print sum}" < out.lt`
      assertTrue "ltemplify.py doesn't like improper type names with spaces" "[ $N_spaces_IMPROPER_NAMES -eq 7 ]"
    cd ../
    rm -rf deleteme/

  cd ../
}

. shunit2/shunit2
