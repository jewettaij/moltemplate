#!/usr/bin/env bash
#
# WHAT THIS SCRIPT DOES:
#
# This script creates the moltemplate (LT) files describing our polymer.
# Then it uses moltemplate.sh to create the LAMMPS input files we need.
#
# PREREQUISITES
#
# The users should have followed the instructions in STEP_1 and created a
# file named "init_crds_polymer_backbone.raw" containing the positions
# of every monomer in the polymer.
#
# PARAMETERS
#
# Make sure that the following numbers are consistent with STEP_1


L_BOND=6.972           # length of the bonds connecting monomers
L_MONOMER=13.944       # length of each monomer (Note: In this model, each
                       # monomer has 3 particles, and two different bonds
                       # along the backbone.  Hence its length equals 2 bonds.
N_TWIST_MOTORS=400     # number of twist motors to insert into the polymer
TORQUE=1.1      # in (kcal/mole)/radian
       # SEE "verbose_version" for an explanation how I estimated this number.



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
  # (However, if you are using twist motors, it should not matter what -helix
  #  parameter you choose, because eventually the supercoiling will be
  #  determined by the amount of torque you apply to these motors.)
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
  #   kB*T=0.001987207*300 (kcal/mole, assuming we are using these energy units)
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
  #  As mentioned, the resulting polymer using these angular potentials
  #  does not have the right mechanical properties (persistence lengths).
  #  After several rounds of iteration, we now use:
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
  #  to prevent the polymer from passing through itself.)

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
  # This was specified in STEP_1 in the $N_MONOMERS variable.  We could ask
  # the user to specify this again here (but they might make a mistake).
  # Instead, we can infer that from the number of non-empty lines in the
  # "init_crds_polymer_backbone.raw" file from STEP_1.  The next line does that:

  N_MONOMERS=`awk '{if ((NF>0) && (substr($1,1,1)!="#")) {n++}} END{print n}' < init_crds_polymer_backbone.raw`

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
    -locations-periodic $N_TWIST_MOTORS 0 \
    -dihedral Disable r c2 c2 r 0 0 1 1 \
-fix-nbody 4 "fix_twist.in" fxTw all twist torque r c2 c2 r 0 0 1 1 "$TORQUE" \
-set-atoms 4 "In Types" "type" r c2 c2 r 0 0 1 1 Rmotor C1motor C1motor Rmotor \
    >> dna_polymer.lt
  # -------------------- END OF TWIST MOTOR SECTION -----------------------



  # The "system.lt" file contains a line that links to the "constraints.lt" file
  # So we must create this file.  By default it is an empty file:
  # (Later on, we will add constraints and run moltemplate again.)
  echo '' > constraints.lt


  # -------------- Finally run moltemplate ------------------

  moltemplate.sh system.lt


  # Now move the files that LAMMPS needs into the directory where we plan
  # to run LAMMPS (the parent directory)

  mv -f system.in* system.data fix_twist*.in delete_link_bonds.in pair_*.in table_*.dat vmd_commands.tcl ../  2> /dev/null
  

  # ---------------- LINKED VERSION -----------------

  # Now create a version of the same circular polymer, but this time
  # add bonds that link the ends of the circular polymer at either end together
  # with the periodic image of the ring polymer on the other side.
  # (In other words, link the origin of replication at one end, with the
  #  terminus of replication at the other end, though the periodic boundary.)

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
  mv -f system.in* system.psf fix_twist*.in delete_link_bonds.in pair_*.in table_*.dat vmd_commands.tcl ../  2> /dev/null
  
  # (We are done building all the files we need for LAMMPS.)

  # -------- OPTIONAL:  Delete temporary files we created earlier ---------

  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate.  They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/

  # Delete the local temporary .lt files we created
  #rm -f init_crds_polymer_backbone.raw
  rm -f dna_polymer.lt dna_forcefield_nb.lt optional_nonbonded.lt constraints.lt

  # Delete any local .pyc files created by running python
  rm -rf __pycache__ *.pyc    2> /dev/null

cd ../
