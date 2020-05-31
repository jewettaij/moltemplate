# PREREQUISITES
#
# The users should have followed the instructions in STEP_1 and created a
# file named "init_crds_polymer_backbone.raw" containing the positions
# of every monomer in the polymer.

# Use these commands to generate the LAMMPS input script and data file

# --- MAKE SURE THESE PARAMETERS AGREE WITH STEP_1 and "dna_forcefield.lt" ---

L_BOND=6.972           # length of the bonds connecting monomers
L_MONOMER=13.944       # length of each monomer (Note: In this model, each
                       # monomer has 3 particles, and two different bonds
                       # along the backbone.  Hence its length equals 2 bonds.


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
  #
  # Notes on the arguments:
  # The "-helix" parameter represents the twist-per-monomer (Δφ) at the
  # start of the simulation.  It is proportional to the "supercoil-density" (σ)
  # The formula for the twist-per-monomer which is:
  #    Δφ = 360*σ*(n_b/10.5)     (in degrees)
  # where n_b is the number of base pairs per monomer, and 10.5 is the
  # natural period of DNA twist in the relaxed state in base-pairs.
  # Example:
  # genpoly_lt.py -helix -72.0 \
  # For details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_lt.md

  genpoly_lt.py \
      -circular yes \
      -helix -14.1545 \
      -bond Backbone c2 c1 \
      -dihedral Backbone r c2 c2 r 0 0 1 1 \
      -polymer-name 'DNAPolymer' \
      -inherits 'DNAForceField'  \
      -monomer-name 'DNAMonomer' \
      -header 'import "dna_monomer.lt"' \
      -padding ${L_BOND},500,500 \
      < init_crds_polymer_backbone.raw > dna_polymer.lt




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

  echo '' >> dna_polymer.lt
  echo 'import "dna_twist_motor.lt"' >> dna_polymer.lt
  echo '' >> dna_polymer.lt

  # How many monomers are there in the polymer?
  # We can infer that from the number of non-empty lines in the
  # "init_crds_polymer_backbone.raw" file.

  N_MONOMERS=`awk '{if ((NF>0) && (substr($1,1,1)!="#")) {n++}} END{print n}' < init_crds_polymer_backbone.raw`

  TORQUE=1.10657  # in (kcal/mole)/radian
                  # SEE BELOW for an explanation how I estimated this number.

  # The "genpoly_modify_lt.py" script will add modifications to an existing
  # polymer created by "genpoly_lt.py", typically at many locations on the
  # polymer (either at evenly spaced intervals or at user specified locations).
  # (Note: This script is distributed and installed along with moltemplate.)
  # Here we use it to add twist motors (and disable any preexisting dihedral
  # interactions that would otherwise prevent those twist motors from spinning).
  # (We also change the atom types for visualization purposes only.)
  # For more details, see:
  # https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_modify_lt.md

  genpoly_modify_lt.py \
    -polymer-name DNAPolymer \
    -length $N_MONOMERS \
    -locations-periodic 400 0 \
    -dihedral Disable r c2 c2 r 0 0 1 1 \
    -fix-nbody 4 "fix_twist.in" fxTw all twist torque r c2 c2 r 0 0 1 1 "$TORQUE" \
    -set-atoms 4 "In Types" "type" r c2 c2 r 0 0 1 1 Rmotor C1motor C1motor Rmotor \
    >> dna_polymer.lt

  # How did we choose the TORQUE (twist motor strength) parameter?
  # This would be the torque necessary to introduce a supercoil
  # density of σ (which is usually measured by experiment to be ~ 0.03-0.05).
  # This is difficult to know because it depends on the shape of the polymer.
  # The torque necessary to introduce a supercoil density σ into a straight
  # polymer made from monomers (length b) connected with torsional springs
  # of strength k_t is
  #  torque = k_t * σ*(2*pi/10.5)*(b/0.332)
  # where σ*(2*pi/10.5) is the Δtwist-per-base-pair corresponding to σ, and
  # (b/0.332) is the number of base pairs per monomer (b in monomer length in nm
  # 0.332nm base pair spacing), and 10.5 is the periodicity of DNA (base pairs)
  # ... AND the spring constant is approximated by:
  #     k_t = (L_t/b)*kB*T
  #     L_t = torsional persistence length of DNA
  #    kB*T = 0.001987207*300kCal/mole,
  # Setting b = 10.5*0.332,  and L_t = 111.0:
  #
  # --> torque = k_t*σ*(2*pi/10.5)*(b/0.332) = 5.875137 (kcal/mole)/ radian
  #     torque * 2*pi / (kB*T) = 61.92 kB*T / turn
  #
  # However, this estimate is for a straight polymer.
  # For a curved polymer that contains both twist and writhe, you can
  # estimate the amount of twist in that polymer by measuring the average 
  # torsional angle between successive monomers in that polymer.
  # First we must generate a curved polymer with a known superhelical density σ
  # and with a conformation that you consider to be realistic.  One way to do
  # this is to create a simple, flat circular polymer (with zero writhe)
  # which is already twisted along its axis by an amount consistent with σ.
  # Then allow the polymer to relax slowly. Hopefully the resulting conformation
  # is similar to the desired conformation.  Then measure the twist stored in
  # the polymer by measuring average torsional angle between monomers.
  # (Alternatively you can measure the polymer's writhe from its shape, and
  #  then subtract writhe from the total supercoil density, if it is known.)
  # From the average torsional angle, you can estimate the torque necessary
  # to create that angle between successive monomers.
  # Then we can insert twist motors into the polymer that exert
  # this amount of torque to achieve the same effect.
  # (This is how I obtained the value of TORQUE used above.
  #  See the file "extract_dihedrals_and_infer_torque.sh" for more details.)
  # -------------------- END OF TWIST MOTOR SECTION -----------------------





  # The "system.lt" file contains a line that links to the "constraints.lt" file
  # So we must create this file.  By default it is an empty file:
  # (Later on, we will add constraints and run moltemplate again.)
  echo '' > constraints.lt





  # Finally run moltemplate

  moltemplate.sh system.lt




  # Now move the files that LAMMPS needs into the directory where we plan
  # to run LAMMPS (the parent directory)

  mv -f system.in* system.data fix_twist*.in delete_link_bonds.in pair_*.in table_*.dat vmd_commands.tcl ../
  
  

  # ---------------- LINKED VERSION -----------------

  # Now create a version of the same circular polymer, but this time
  # add bonds that link the ends of the circular polymer at either end together
  # with the periodic image of the ring polymer on the other side.
  # (In other words, link the origin of replication at one end, with the
  #  terminus of replication at the other end, though the periodic boundary.)
  # This is one of several strategies one could use to keep the polymer under
  # tension as we let the system collapse (contract) using "fix deform".
  # We don't want to just completely let go of the polymer ends, because
  # the polymer needs to be under some kind of tension for the plectonemic
  # supercoils to slowly form along it's length correctly as the tension is
  # relaxed.  Furthermore the polymer needs to be linked in 2 places (not 1)
  # to the polymers on the other side of the periodic boundary.  This is to 
  # prevent the ends of the circle from spinning freely.  Otherwise, the
  # entire circular polymer would simply twist around its own long axis (the
  # x-axis) and form one incredibly long plectonemic supercoil.  (This is not
  # what bacterial nucleoids look like.)  Linking the polymer in 2 nearby
  # places along the chain avoids this.  And this also does a more accurate job
  # of simulating the replication of DNA in Caulobacter crescentus, because
  # in that bacteria either end of the cicular DNA (ie the Ori and Ter regions)
  # are anchored to opposite sides of the cell wall (in a way that would prevent
  # them from freely rotating).  Since there is no cell wall in this model,
  # the only thing left to bind the polymer to its periodic image.

  N_MONOMERS_HALF=`awk -v N=$N_MONOMERS 'BEGIN{print int(N/2)}'`
  N1=0
  N2=$(( N_MONOMERS_HALF-1 ))
  N3=$N_MONOMERS_HALF
  N4=$(( N_MONOMERS-1 ))

  echo 'DNAPolymer {' > constraints.lt
  echo '  #(augment the definition of "DNAPolymer" by adding 2 bonds)' >> constraints.lt
  echo '  write("Data Bonds") {' >> constraints.lt
  echo "    \$bond:bond_link1 @bond:Periodic \$atom:mon[$N1]/c1 \$atom:mon[$N2]/c1" >> constraints.lt
  echo "    \$bond:bond_link2 @bond:Periodic \$atom:mon[$N3]/c1 \$atom:mon[$N4]/c1" >> constraints.lt
  echo '  }' >> constraints.lt
  echo '' >> constraints.lt
  echo '  write_once("In Settings") {' >> constraints.lt
  echo "    # type of bond used to connect 2 ends of the polymer together" >> constraints.lt
  echo "    bond_coeff @bond:Periodic harmonic 10.0 $L_BOND" >> constraints.lt
  echo '  }' >> constraints.lt
  echo '}' >> constraints.lt





  # Now run moltemplate again with these new constraints added

  moltemplate.sh system.lt




  # Now move the files that LAMMPS needs into the directory where we plan
  # to run LAMMPS (the parent directory), but change its name so that
  # we don't overwrite the "system.data" file we created the last time.
  # (Later, we will need both these files.)

  mv -f system.data ../system_linked.data

  # Move the other files to the parent directory as well (but don't rename them)
  mv -f system.in* system.psf fix_twist*.in delete_link_bonds.in pair_*.in table_*.dat vmd_commands.tcl ../
  



  # (We are done building all the files we need for LAMMPS.)






  # -------- OPTIONAL:  Delete temporary files we created earlier ---------

  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate.  They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/

  # Delete the local temporary .lt files we created
  rm -f dna_polymer.lt dna_forcefield_nb.lt optional_nonbonded.lt init_crds_polymer_backbone.raw

  # Delete any local .pyc files created by running python
  rm -rf __pycache__ *.pyc

  # If we created a temporary python environment, remove that now:
  #deactivate
  #rm -rf venv

cd ../










  # ---------------------------- PLEASE IGNORE --------------------------------
  # ---------------------------- COMMENTED OUT --------------------------------
  # ---- OLD METHOD (does not work well for beads larger than Debye length) ---
  ## Create the tabular pair tables that the polymer needs:
  ##Define interactions between the central particles that make up the DNA chain
  #
  ##echo "" > dna_forcefield_nb.lt
  #./gen_nonbonded_tables.py \
  #    -types "DNAForceField/C* DNAForceField/C*" \
  #    -label "DNA_U0=inf" \
  #    -table-file "table_dna_U0=inf.dat" \
  #    -Ntable 64 \
  #    -script-file "In Settings" \
  #    -pair-style "table" \
  #    -barrier-height inf \
  #    -lj 1.0 6.972 \
  #    -lambda 1 \
  #    -ljrmax 6.972 \
  #    -charge -42.0 \
  #    -Ldebye 0.000000001 \
  #    -Ldebye-cut 0.0 \
  #    -escale 0.5961621 \
  #    >> dna_forcefield_nb.lt
  #
  # Incidentally, "0.5961621" = kB*T in kCal/mole (assuming T=300K)
  #
  ## Define the interaction between the DNAForceField/R particles with all
  ## other particles ("*"). Set the barrier-height to 0 to turn off these forces
  #
  #echo "" >> dna_forcefield_nb.lt
  #./gen_nonbonded_tables.py \
  #    -types "DNAForceField/R *" \
  #    -label "DNA_U0=0" \
  #    -table-file "table_dna_U0=0.dat" \
  #    -Ntable 64 \
  #    -script-file "In Settings" \
  #    -pair-style "table" \
  #    -barrier-height 0 \
  #    -lj 1.0 6.972 \
  #    -lambda 1 \
  #    -ljrmax 6.972 \
  #    -charge -42.0 \
  #    -Ldebye 0.000000001 \
  #    -Ldebye-cut 0.0 \
  #    -escale 0.5961621 \
  #    >> dna_forcefield_nb.lt
  #
  ##echo "} # Polymer (pair tables)" >> dna_polymer.lt
  #
  ## Now generate the files defining the OPTIONAL forces between particles
  ## used in this DNA model.  We do this for a range of parameters
  ## that we might use (a range of U0 barrier height values).
  ## Later we can select from these different models by selecting a different
  ## script file (beginning with "pair_...")
  #rm -f optional_nonbonded.lt
  #for BARRIER_HEIGHT in inf 64 32 16 8 4 2 1 0.5 0; do
  #  ./gen_nonbonded_tables.py \
  #     -types "DNAForceField/C* DNAForceField/C*" \
  #     -label "DNA_U0=${BARRIER_HEIGHT}kT" \
  #     -table-file "table_dna_U0=${BARRIER_HEIGHT}kT.dat" \
  #     -Ntable 64 \
  #     -script-file "pair_dna_U0=${BARRIER_HEIGHT}kT.in" \
  #     -pair-style "table" \
  #     -barrier-height $BARRIER_HEIGHT \
  #     -lj 1.0 6.972 \
  #     -lambda 1 \
  #     -ljrmax 6.972 \
  #     -charge -42.0 \
  #     -Ldebye 0.000000001 \
  #     -Ldebye-cut 0.0 \
  #     -escale 0.5961621 \
  #     >> optional_nonbonded.lt
  #  # ("inf"<-->infinity, which is implemented here as a large positive number)
  #done
  #
  ## Notes:
  ##
  ##  By default, the energy between particles is defined as:
  ##      U(r) = U_LR(r) + U_yuk(r)
  ## -lj supplies the Lennard-Jones parameters (ε, σ)
  ##
  ## -lambda (λ) supplies the (optional) 3rd Lennard-Jones parameter.
  ##   --> U_LR(r) = ε * ((σ/(r-rshiftLJ))^12 - 2*λ*(σ/(r-rshiftLJ))^6)
  ## -charge = the charge of the particle
  ##
  ## -Ldebye = the Debye length (electrostatic decay length due to counter ions)
  ##   --> U_yuk(r) = (Ke/r) * exp(-(r-rshiftQ)/Ldebye)     ("yukawa potential")
  ##             Ke = ke*qi*qj / eps_r
  ##             ke = 8.9875517873681764e09 J*m*C^-2 = 1/4*pi*eps_0
  ##                  (https://en.wikipedia.org/wiki/Coulomb's_law)
  ##
  ## -barrier-height U0   If this parameter is specified, then U(r) changes to
  ##            U(r) = (U0/(pi/2)) * arctan(Uorig(r)/U0)
  ##  where Uorig(r) = the original U(r) potential.  
  ##                  This function approaches U0 as r->0 (a non-infinite value)
  ##
  ## -escale multiplies the energy parameters (ε,U0) by a constant before use.
  ##         (Typically I choose that constant to be kB*Temperature)
  ##
  ## - Why did I choose Ldebey ≈ 0 and Ldebye-cut = 0? -
  ## Because the range of the Lennard-Jones spheres we are using to prevent the
  ## polymer from passing through itself is wider than the Debye length there
  ## is no longer a reason to consider electrostatic repulsion in this model.
  ## Doing so would just make the width of the polymer even larger.  So in this
  ## case, I turned off the electrostatics by setting Ldebye to a tiny number,
  ## and Ldebye-cut to 0.  (Note: For typical salt concentrations Ldebye ≈ 1nm.)
  #
  # ---------------------------- PLEASE IGNORE --------------------------------
