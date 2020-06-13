# Use these commands to generate the LAMMPS input script and data file


# First, make sure that the following constants are consistent with STEP_1
L_MONOMER=0.98293    # relaxed distance between monomers along the polymer
N_IHF=1070           # number of IHF proteins bound to the polymer
N_TWIST_MOTORS=40    # number of twist motors distributed along the polymer
WIDTH_IHF=10         # how many monomers is IHF bound to?
WIDTH_TWIST_MOTORS=4 # how many monomers are occupied by each twist motor?
TORQUE=0.4

# N_MONOMERS = the number of monomers in the polymer
#            = the number of lines in the "init_crds_polymer_backbone.raw" file
N_MONOMERS=`awk '{if ((NF>0) && (substr($1,1,1)!="#")) {n++}} END{print n}' < moltemplate_files/init_crds_polymer_backbone.raw`

# ---- Optional: Shorten the physical length of the polymer. ----
#           We do this because the polymer has N_IHF bend modifications, and
#           each one reduces the polymer length by roughly 2*WIDTH_IHF monomers.
#           Reducing the length of L_MONOMER will cause the monomers to be
#           closer together.  (We will minimize the system before the main
#           simulation to prevent steric collisions).  I think reducing the
#           overall polymer length will improve numerical stability, by not
#           overstretching the polymer at the beginning of the simulation.
#           (But I don't think it is strictly necessary to do this.)

L_MONOMER=`awk -v N=$N_MONOMERS -v L=$L_MONOMER -v n=$N_IHF -v w=$WIDTH_IHF 'BEGIN{print L*(N-n*w*2)/N}'`

# ---- Optional ----



# Now, create LAMMPS input files this way:

cd moltemplate_files

  # Use the "genpoly_lt.py" to generate a moltemplate file (.LT file)
  # describing the polymer you want to simulate.  You must specify the
  # name of the moltemplate object which will be used as the monomer subunit
  # in the final polymer (eg. "DNAMonomer"), as well as any bonds (or angles
  # or dihedrals) linking one monomer to the next monomer, as well as the
  # helical twist angle (if applicable).  All of the details regarding
  # the behaviour of the polymer are contained in the "dnamonomer.lt" file
  # which defines the "DNAMonomer" object, as well as a link to the file
  # which defines "DNAForceField" (which DNAMonomer uses).  For details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_lt.md

  genpoly_lt.py -circular yes \
                -helix 102.7797 \
                -bond Backbone a a \
                -bond Backbone b b \
                -dihedral MajorGroove b b a a 0 1 1 2 \
                -dihedral Torsion a a b b 1 0 0 1 \
                -polymer-name 'DNAPolymer' \
                -inherits 'DNAForceField'  \
                -monomer-name 'DNAMonomer' \
                -header 'import "dna_monomer.lt"' \
                -padding 2,100,100 \
                < init_crds_polymer_backbone.raw > dna_polymer.lt

  # (Note: The "-helix" parameter represents the twist-per-monomer (Δφ) at the
  #        start of the simulation.  Example "genpoly_lt.py -helix 102.857 ...")

  
  # This polymer is bound to hundreds of IHF proteins and twist motors.  Use
  # "import" to load the definition of the "IHFProtein" molecule we will add.

  echo '' >> dna_polymer.lt
  echo 'import "ihf_protein.lt"' >> dna_polymer.lt
  echo '' >> dna_polymer.lt

  # We have to distribute both the IHFProtein and twist motors on the DNA in a
  # random way so that they don't overlap with each other. That is what the next
  # two commands do.  The "genpoly_modify_lt.py" program (which is distributed
  # with moltemplate) can be used to generate lists of integers which avoid
  # overlap with each other.  If, for some reason, you want more details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_modify_lt.md

  # Generate a file containing the locations of the twist motors on the polymer
  genpoly_modify_lt.py \
      -circular yes \
      -length $N_MONOMERS \
      -locations-periodic $N_TWIST_MOTORS 0 \
      -width $WIDTH_TWIST_MOTORS \
      -write-occupancy dna_polymer_occupancy.txt \
      -write-locations locations_twist_motors.txt

  # Generate a file containing the locations of the IHF proteins on the polymer
  # (taking care not to overlap with the twist motors)
  genpoly_modify_lt.py \
      -circular yes \
      -length $N_MONOMERS \
      -locations-random $N_IHF 0 \
      -locations-random-attempts 50 \
      -width $WIDTH_IHF \
      -read-occupancy dna_polymer_occupancy.txt \
      -write-locations locations_proteins.txt

  # NOTE: We haven't actually modified the polymer yet.
  #       All we did was decide where we want the modifications to go.
  #       In the commands below, we will modify the polymer.


  # -------- modify the polymer, adding the IHF proteins ----------
  # Loop over the integers in the "locations_proteins.txt" file.
  # Each line of that file contains an integer (=$i) indicating (approximately)
  # the monomer where the IHF protein is bound


  # These modifications are to the "DNAPolymer" object, so enclose them
  # in "DNAPolymer {" ... "}"
  echo "" >> dna_polymer.lt
  echo "DNAPolymer {" >> dna_polymer.lt
  echo "  #(augment the definition of the \"DNAPolymer\" object)">>dna_polymer.lt
  echo "" >> dna_polymer.lt


  I=0
  for i in $(cat locations_proteins.txt); do
    # Variables:
    #
    # I is an integer from 0 to the number of IHF proteins-1
    # which represents the IHF protein you are considering.
    #
    # i is an integer (from 0 to the number of $N_MONOMERS-1)
    # which indicates the first monomer from the polymer that
    # this IHF protein is bound to.
    #
    # The position of other monomers in the polymer is stored in
    # ip1, ip4, ip5, ip8, which are usually i+1, i+4, i+5, i+8.
    # However since i ranges from 1 to $N_MONOMERS, it is possible
    # that these numbers could exceed $N_MONOMERS.  Since this is
    # a circular polymer, we take the remainder of i+1 (and i+4,...)
    # after division by $N_MONOMERS using the % operator.
    ip1=$(( i+1 % $N_MONOMERS ))
    ip4=$(( i+4 % $N_MONOMERS ))
    ip5=$(( i+5 % $N_MONOMERS ))
    ip8=$(( i+8 % $N_MONOMERS ))
    ip9=$(( i+9 % $N_MONOMERS ))
    # Add a bond directly between two atoms in the DNA
    # This bond is what causes the "kink" in the shape of the DNA.
    echo "" >> dna_polymer.lt
    echo "  write(\"Data Bonds\") {" >> dna_polymer.lt
    echo "    \$bond:IHFDD_${I} @bond:IHFProtein/DNA_DNA_IHF \$atom:mon[$i]/a \$atom:mon[$ip9]/b" >> dna_polymer.lt
    echo "  }" >> dna_polymer.lt
    # OPTIONAL: Create a new particle and attach it to several
    #           atoms in this part of the polymer.  The purpose of
    #           that particle is for the purpose of visualization.
    #           We want to display a particle at the location of
    #           where the actual IHF protein would have been.
    echo "" >> dna_polymer.lt

    # Create the new particle and move it to somewhere (approximately)
    # near where the corresponding monomer should be that it will bind to.
    # (It doesn't matter if this particle is not in the optimimal place because
    #  we will minimize the system before simulating it.  Minimization will
    #  relax any stretched bonds and reduce steric overlap that might occur.)

    # The x position should be (i+WIDTH_IHF/2)*L_MONOMER
    # (I also subtract by the half the length of the polymer to center it.)
    XPOS=`awk  -v i=$i  -v L=$L_MONOMER  -v w=$WIDTH_IHF  -v n=$N_MONOMERS  'BEGIN{print (i+w/2-n/2)*L}'`
    # The y and z positions should be close (but not too close)
    # to the polymer axis.
    YPOS=`awk -v L=$L_MONOMER 'BEGIN{print 4*L}'` #(move it to avoid overlap)
    ZPOS=0
    echo "  ihf[$I] = new IHFProtein.move($XPOS,$YPOS,$ZPOS)" >> dna_polymer.lt
    echo "" >> dna_polymer.lt

    # Now attach the new particle to the polymer with bonds

    echo "  write(\"Data Bonds\") {" >> dna_polymer.lt
    echo "    \$bond:IHF_DNA_${I}_1 @bond:IHFProtein/IHF_DNA \$atom:ihf[$I]/x \$atom:mon[$ip1]/b" >> dna_polymer.lt
    echo "    \$bond:IHF_DNA_${I}_4 @bond:IHFProtein/IHF_DNA \$atom:ihf[$I]/x \$atom:mon[$ip4]/a" >> dna_polymer.lt
    echo "    \$bond:IHF_DNA_${I}_5 @bond:IHFProtein/IHF_DNA \$atom:ihf[$I]/x \$atom:mon[$ip5]/b" >> dna_polymer.lt
    echo "    \$bond:IHF_DNA_${I}_8 @bond:IHFProtein/IHF_DNA \$atom:ihf[$I]/x \$atom:mon[$ip8]/a" >> dna_polymer.lt
    echo "  }" >> dna_polymer.lt

    I=$(( I+1 ))
  done

  echo "} #DNAPolymer" >> dna_polymer.lt
  echo "" >> dna_polymer.lt



  # Note: if you don't care about adding extra particles for visual reasons
  #       you can try replacing the entire for loop above with the following:
  #       This will place a single bond between mon[i]/a and mon[i+9]/b
  #
  # genpoly_modify_lt.py \
  #     -circular yes \
  #     -polymer-name DNAPolymer \
  #     -length $N_MONOMERS \
  #     -locations locations_proteins.txt \
  #     -bond IHF a b 0 9 \
  #     >> dna_polymer.lt



  # Add twist motors.
  # If I only wanted to add a single twist motor, it would be easy to manually
  # add some extra lines to the "dna_polymer.lt" file.  However here I wrote
  # this script to make it possible to put many, many twist motors along the
  # polymer.  To do that, I created a new script named "genpoly_modify_lt.py"
  # which generates many modifications to a polymer at user-defined locations.
  # It's overkill for what we need in this example since we only use 1 motor.

  # "genpoly_modify_lt.py" needs to know the length of the polymer we created.
  # Count the number of non-blank, non-comment lines in the coordinate file:

  echo '' >> dna_polymer.lt
  echo 'import "dna_twist_motor.lt"' >> dna_polymer.lt
  echo '' >> dna_polymer.lt

  genpoly_modify_lt.py \
      -circular yes \
      -polymer-name DNAPolymer \
      -length $N_MONOMERS \
      -locations locations_twist_motors.txt \
      -bond Motor   a a 1 2 \
      -bond Disable b b 1 2 \
      -dihedral MajorGrooveML b b a a 0 1 1 2 \
      -dihedral MajorGrooveMR a a b b 1 2 2 3 \
      -dihedral Disable       a a b b 2 1 1 2 \
      -dihedral Disable       b b a a 1 2 2 3 \
      -dihedral Disable       b a a b 1 1 2 2 \
      -set-atoms 6 "system.in.types" "type" b b a a b b 0 1 1 2 2 3 Bm2 Bm Am Am Bm Bm2 \
      -fix-nbody 4 "fix_twist.in" fxTw all twist torque b a a b 1 1 2 2 "$TORQUE" \
      >> dna_polymer.lt

  # NOTE: To force the motor to twist at a constant rate (instead of applying
  # a constant torque), use this instead.
  #
  # -fix-nbody 4 "fix_twist_rate_5.0_100_14400.in" fxTw all twist torque b a a b 1 1 2 2 "5.0 100 14400"
  # (WARNING:  Simulation can become numerically unstable if twisted too far.)




  # ---------- OPTIONAL -------------------------------------
  # --- Delete the bond interfering with the twist motor. ---
  echo '' >> dna_polymer.lt
  echo 'DNAPolymer {' >> dna_polymer.lt
  # Note: We already disabled this bond using "-bond Disable b b 1 2"
  #       (by setting its spring constant to 0).  However you actually have
  #       to delete that bond if you want it not to appear in visualization
  #       software tools like VMD (which was my goal).  To delete the bond,
  #       you have to know its $bond: name.  Bonds generated by genpoly_lt.py
  #       have names like "genp_bondi_j", where "j" indicates the monomer (from
  #       mod_locations.txt) and "i" represents the bond-per-monomer (2 here).
  awk -v N=$N_MONOMERS '{print "  delete genp_bond2_"1+($1+1)%N}' < locations_twist_motors.txt >> dna_polymer.lt
  awk -v N=$N_MONOMERS '{print "  delete gpm_bond2_"1+($1)%N}' < locations_twist_motors.txt >> dna_polymer.lt
  echo '}' >> dna_polymer.lt
  # ---------- OPTIONAL -------------------------------------

  # Now clean up by discarding temporary files we no longer need:
  rm -f dna_polymer_occupancy.txt locations_*.txt




  # ------------- Run moltemplate.sh -------------

  # Then run moltemplate.sh on "system.lt".
  # (Note: "system.lt" contains a reference to the polymer file we created.)

  moltemplate.sh system.lt

  # This will generate various files with names ending in *.in* and *.data. 
  # These files are the input files directly read by LAMMPS.  Move them to 
  # the parent directory (or wherever you plan to run the simulation).
  mv -f system.in* fix_twist*.in system.data ../

  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate. They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/

  # Optional: Delete other temporary files:
  rm -f init_crds_polymer_backbone.raw
  rm -f dna_polymer.lt

cd ../

