# Use these commands to generate the LAMMPS input script and data file


# Create LAMMPS input files this way:

cd moltemplate_files

  # Use "genpoly_lt.py" to create a polymer which follows the curve 
  # traced by the coordinates that we copied into this directory earlier.
  # (Note: "genpoly_lt.py" is distributed and installed along with moltemplate.)
  # The "genpoly_lt.py" program generates a moltemplate file (.LT file)
  # describing the polymer you want to simulate.  You must specify the
  # name of the moltemplate object which will be used as the monomer subunit
  # in the final polymer (eg. "DNAMonomer"), as well as any bonds (or angles
  # or dihedrals) linking one monomer to the next monomer, as well as the
  # helical twist angle (if applicable).  All of the details regarding
  # the behaviour of the polymer are contained in the "dnamonomer.lt" file
  # which defines the "DNAMonomer" object, as well as a link to the file
  # which defines "DNAForceField" (which DNAMonomer uses).

  # Notes on the arguments:
  # The "-helix" parameter represents the twist-per-monomer (Δφ) at the
  # start of the simulation.  It is proportional to the "supercoil-density" (σ)
  # The formula for the twist-per-monomer which is:
  #    Δφ = 360*σ*(n_b/10.5)     (in degrees)
  # where n_b is the number of base pairs per monomer, and 10.5 is the
  # natural period of DNA twist in the relaxed state in base-pairs.
  # Example:
  # genpoly_lt.py -helix -14.1545 \
  # For details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_lt.md

  genpoly_lt.py \
      -circular yes \
      -helix 0 \
      -bond Backbone c2 c1 \
      -dihedral Backbone r c2 c2 r 0 0 1 1 \
      -polymer-name 'DNAPolymer' \
      -inherits 'DNAForceField'  \
      -monomer-name 'DNAMonomer' \
      -header 'import "dna_monomer.lt"' \
      -padding 30,30,110 \
      < init_crds_polymer_backbone.raw > dna_polymer.lt

  # (Note: The "-helix" parameter represents the twist-per-monomer (Δφ) at the
  #        start of the simulation.  Example "genpoly_lt.py -helix 102.857 ...")




  # --------- angles and dihedrals -------------
  #
  # We are using angle and dihedral interactions who's energy (near the minima)
  # is approximated by:
  #
  #    Ubend(theta) = (k_b/2)*theta^2   
  #    Utorsion(phi) = (k_t/2)*phi^2   
  #         (k_b and k_t in energy/rad^2)
  #
  # For small b, k_b and k_t are proportional to the bending and torsional
  # persistence lengths (L_b and L_t) according to:
  #
  #   L_b/b ~= k_b / kB*T
  #   L_t/b ~= k_t / kB*T
  #
  #   b = n_b * 0.332nm
  #   n_b=the number of base pairs per monomer (typically less than 50)
  #   kB*T=0.001987207*300 (kCal/mole, assuming we are using these energy units)
  #
  # However for numerical stability we must forbid certain angles from being
  # visited, and this requires changing the potential to something slightly more
  # complicated.  The new Ubend(angle) and Utorsion(phi), now have this form:
  #
  # U(theta) = ((2/pi)t_range)^2 * K (-1+1/ cos((theta-theta0)/((2/pi)t_range)))
  #
  # where t_range is the range of allowable angles around the minima at theta0.
  # (See calc_table_angle.py for details.)
  #
  # As a result, the bending and torsional persistence lengths
  # are no longer exactly equal to k_b/kB*T and k_t/kB*T (respectively).
  # To achieve the desired persistence lengths, we must try several different
  # k_b and k_t values and iterate until the desired persistence lengths are
  # reproduced.  The initial iteration (iteration 0) uses the estimate above.
  #
  #  -------- ITERATION 0: ------
  # 
  #   L_b = 50nm         ->   k_b = (L_b/b)*kB*T = 2.1059845
  #   L_t = 100nm        ->   k_t = (L_t/b)*kB*T = 4.211969
  #   (NOTE: In iteration 0, I set L_t=100.0 by mistake.  Later, in ITERATION1
  #    after reading Gore++Bustamante_Nature_2006, I set L_t = 111nm)
  #
  #./calc_table_angle.py 180.0 2.1059845 120.0 181 CCC 0.0 180.0 65.0 12.0\
  #                      "EQ 0.0"  > table_angle.dat
  #./calc_table_angle.py 0 4.211969 180.0 180 RCCR -179 179 -160 12.0 160 -12.0\
  #                      "DEGREES" > table_dihedral.dat
  #
  #  -------- ITERATION 1: ------
  #
  # L_b_target = 35.0,     L_b_sim1 =  73.77234648
  # L_t_target = 111.0,    L_t_sim1 =  108.046530763
  #
  # k_b_sim1 = k_b_sim0 * (L_b_target / L_b_sim0) = 0.9991475
  # k_t_sim1 = k_t_sim0 * (L_t_target / L_t_sim0) = 4.327104
  #
  #./calc_table_angle.py 180.0 0.9991475 120.0 181 CCC 0.0 180.0 65.0 12.0\
  #                      "EQ 0.0"  > table_angle.dat
  #./calc_table_angle.py 0 4.327104 180.0 180 RCCR -179 179 -160 12.0 160 -12.0\
  #                      "DEGREES" > table_dihedral.dat
  #
  #                        :
  #                        :
  #  -------- SKIPPING ITERATIONS 2 THROUGH 4 ---------
  #                        :
  #                        :
  #
  #  -------- ITERATION 5: ------
  #
  # L_b_target = 50.0,     L_b_sim4 = 52.3242566
  # L_t_target = 111.0,    L_t_sim4 = 104.647553
  #
  # k_b_sim5 = k_b_sim4 * (L_b_target / L_b_sim4) = 1.14983
  # k_t_sim5 = k_t_sim4 * (L_t_target / L_t_sim4) = 4.42169
  #
  # (In iteration 5, I also added a hard-sphere steric (Lennard-Jones)
  #  repulsion to the backbone DNA particles at short distances
  #  to prevent the polymer from passing through itself.  Earlier
  #  iterations did not use this approach.)

  ./calc_table_angle.py 180.0 1.14983 120.0 181 CCC 0.0 180.0 65.0 12.0\
                        "EQ 0.0"  > table_angle.dat
  ./calc_table_angle.py 0 4.42169 180.0 180 RCCR -179 179 -160 12.0 160 -12.0\
                        "DEGREES" > table_dihedral.dat

  # --> L_b_sim5 = 48.6223
  # --> L_t_sim5 = 109.9063




  # Finally add a new table to "table_dihedrals.dat" with the forces turned off:
  ./calc_table_angle.py 0.000 0.000000 180.2 180 ZEROS -179 179 \
                        >> table_dihedral.dat




  # ------------------------- TWIST MOTORS ---------------------------  
  # If you don't want twist motors added to the polymer,
  # comment out this entire section.

  # How many monomers are there in the polymer?
  # We can infer that from the number of non-empty lines in the
  # "init_crds_polymer_backbone.raw" file.

  N_MONOMERS=`awk '{if ((NF>0) && (substr($1,1,1)!="#")) {n++}} END{print n}' < init_crds_polymer_backbone.raw`

  TORQUE=1.10657 # in (kcal/mole)/radian

  # If I only wanted to add a single twist motor, it would be easy to manually
  # add some extra lines to the "dna_polymer.lt" file.  However here I wrote
  # this script to make it possible to put many, many twist motors along the
  # polymer.  To do that, I created a new script named "genpoly_modify_lt.py"
  # which generates many modifications to a polymer at user-defined locations.
  # It's overkill for what we need in this example since we only use 1 motor.

  echo '' >> dna_polymer.lt
  echo 'import "dna_twist_motor.lt"' >> dna_polymer.lt
  echo '' >> dna_polymer.lt

  # Now run the script that makes (potentially)
  # many modifications to the polymer.
  # In our case it will modify the polymer to add a twist motor.
  # The position of that motor is in the file "mod_locations.txt"
  # (which currently only has one entry).  For more details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_modify_lt.md

  genpoly_modify_lt.py \
    -polymer-name DNAPolymer \
    -length $N_MONOMERS \
    -locations mod_locations.txt \
    -dihedral Disable r c2 c2 r 0 0 1 1 \
    -fix-nbody 4 "fix_twist.in" fxTw all twist torque r c2 c2 r 0 0 1 1 "$TORQUE" \
    -set-atoms 4 "In Types" "type" r c2 c2 r 0 0 1 1 Rmotor C1motor C1motor Rmotor \
    >> dna_polymer.lt

  # NOTE: To force the motor to twist at a constant rate (instead of applying
  # a constant torque), use this instead.
  #
  # -fix-nbody 4 "fix_twist_rate_5.0_100_14400.in" fxTw all twist torque b a a b 1 1 2 2 "5.0 100 14400"
  # (WARNING:  Simulation can become numerically unstable if twisted too far.)
  # -------------------- END OF TWIST MOTOR SECTION -----------------------


  # Then run moltemplate on "system.lt".
  # (Note: "system.lt" contains a reference to the polymer file we created.)

  moltemplate.sh system.lt

  # This will generate various files with names ending in *.in* and *.data. 
  # These files are the input files directly read by LAMMPS.  Move them to 
  # the parent directory (or wherever you plan to run the simulation).
  mv -f system.in* fix_twist*.in system.data table_*.dat vmd_commands.tcl ../

  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate. They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/

  # Optional: Delete other temporary files:
  rm -f init_crds_polymer_backbone.raw
  rm -f dna_polymer.lt

cd ../

