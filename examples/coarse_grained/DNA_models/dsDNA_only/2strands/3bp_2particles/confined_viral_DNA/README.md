confined DNA example
==========

### Images

<img src="http://moltemplate.org/images/DNA/3bp2p/3bp2p_dna_monomer_LR.jpg" width=220> <img src="images/plus.svg" height=80> <img src="http://moltemplate.org/images/DNA/3bp2p/HIV_capsid+DNA/curve_ndmansfield_11x11x11_white_LR.jpg" width=200> <img src="images/rightarrow.svg" height=80> <img src="http://moltemplate.org/images/DNA/3bp2p/HIV_capsid+DNA/dna_t=0_green+cyan_scale0.5_width0.2_bbk_occ_light2_LR.jpg" width=200> <img src="images/plus.svg" height=80>
<img src="http://moltemplate.org/images/DNA/3bp2p/HIV_capsid+DNA/capsid_bbk_occ_light3_LR.jpg" width=200> <img src="images/rightarrow.svg" height=80>
<img src="http://moltemplate.org/images/DNA/3bp2p/HIV_capsid+DNA/dna+capsid_t=0_scale0.5_green+cyan_bbk_occ_light2_LR.jpg" width=200> <img src="images/rightarrow.svg" height=80> **LAMMPS**  <img src="images/rightarrow.svg" height=80>  <img src="http://moltemplate.org/images/DNA/3bp2p/HIV_capsid+DNA/dna+capsid_t=8680000_green+cyan_bbk_occ_light2_LR.jpg" width=200>


## Description

This is an example demonstrating how to confine a polymer (eg. DNA) within an
arbitrary shaped box (eg. the mature HIV viral capsid, shown in grey).
Initially the polymer is confined to a rectangular box
which is small enough to fit within the desired container (the viral capsid).
The final two images show the initial and final conformation of
the DNA polymer (before and after the relaxation simulation).
This example uses the ["3bp2p"](../simple_dna_example)
coarse-grained DNA model.


## Overview

Initially, the conformation of the polymer is created with a combination
of [ndmansfield](https://github.com/jewettaij/ndmansfield)
(a random self-avoiding-curve generator), smoothed with
[interpolate_curve.py](../../../../../../../doc/doc_interpolate_curve.md),
and then re-scaled so that the spacing between points along the
curve matches the desired spacing between monomers in the polymer.
(These commands can be entered into a terminal running BASH.)

```shell
  NMONOMERS=3058;           # number of monomers
  B=0.996                   # physical length per monomer
  Nx=11; Ny=11; Nz=11       # number of lattice sites

  # Generate a lattice polymer shape (integer coordinates)
  #   Note: You can download the "ndmansfield" source code here:
  #   https://github.com/jewettaij/ndmansfield
  ./ndmansfield -box $Nx $Ny $Nz -seed 0 -tstart 1 -tsave 1000000 -tstop 1000000 > curve_lattice.txt

  # Interpolate and rescale the coordinates:
  SCALE=`echo "($NMONOMERS*$B)/($Nx*$Ny*$Nz-1)" | bc -l`   #physical curve length / #lattice sites
  interpolate_curve.py $NMONOMERS $SCALE < curve_lattice.txt > curve_smooth.txt
```

Then the
[genpoly_lt.py](../../../../../../../doc/doc_genpoly_lt.md),
tool is used to create a moltemplate file ("dna_polymer.lt")
which defines a polymer which lies along the curve
created in step 1 (see above).

```shell
  # Create a moltemplate file ("dna_polymer.lt") for a polymer with this shape
  genpoly_lt.py -helix 102.7797 \
                -polymer-name 'DNAPolymer' \
                -monomer-name 'DNAMonomer' \
                -inherits 'DNAForceField' \
                -header 'import "dna_monomer.lt"' \
                -bond Backbone a a \
                -bond Backbone b b \
                -dihedral MajorGroove b b a a 0 1 1 2 \
                -dihedral Torsion a a b b 1 0 0 1 \
                -padding 20,20,20 \
                < curve_smooth.txt > dna_polymer.lt
```

Finally moltemplate is invoked to create a LAMMPS data file and input script
containing the polymer (as well as the box that it is contained within).

```shell
  moltemplate.sh system.lt
```

(Note: The "system.lt" file contains a link to "dna_polymer.lt".)

For a detailed step-by-step explanation of this process
(with an abundance of crufty comments),
see these two files:
["STEP_1_generate_coords.sh"](STEP_1_generate_coords.sh)
and
["STEP_2_generate_LAMMPS_files.sh"](STEP_2_generate_LAMMPS_files.sh).


Some of the particles in the simulation form the walls of the container
(the viral capsid shell).  The are immobilized by excluding them from
the group of atoms that is supplied to the integrator.  (For details, see
"run.in.min", and pay attention to the "group" and "fix nve/limit" commands.)


##    Prerequisites

LAMMPS must be compiled with the "MOLECULE" AND "USER-MISC" packages enabled.
If you receive this error message (or something similar):
"dihedral_style spherical: Unknown dihedral style", then you must follow
[these instructions](https://lammps.sandia.gov/doc/Build_package.html),
and recompile LAMMPS.


## Details

This example uses the "3bp2p" DNA model described [here](../simple_dna_example).
The force field parameters were tuned to reproduce realistic geometry and
mechanical properties of DNA, including a major/minor groove, helicity, length,
persistence length, and torsional persistence length.


## Polymer Melts

Note: In this example, there was only 1 polymer, but you can create a
system with multiple polymers (of various lengths) confined in the same box
by running the "genpoly_lt.py" script with the "-cuts" argument in step 2.
This is a useful way to generate polymer melts.  However if the polymers
in your simulation are long, then their conformation is likely to be
kinetically trapped due to entanglement.  In that case, you must use
one of the following methods below to relax the polymer to its
equilibrium conformation.


## Warning: This is not a realistic model of DNA in the HIV virus

### Soft potentials needed to escape kinetic traps

It is impossible to equilibrate a long confined polymer using this method.
The duration of this simulation is insufficient for a polymer
of this length to reach a realistic equilibrium conformation.
In this example, the DNA is expected to have a liquid-crystal-like
equilibrium conformation, with most of the DNA polymers aligned
parallel to each other.  But in this short simulation, the DNA
remains bent and entangled with itself,

In order to reach a plausible equilibrium conformation,
it is necessary to allow the polymer to pass through itself.
In LAMMPS the later can be done using any of the following pair styles:

https://lammps.sandia.gov/doc/pair_fep_soft.html

https://lammps.sandia.gov/doc/pair_table.html

https://lammps.sandia.gov/doc/pair_soft.html

https://lammps.sandia.gov/doc/pair_gauss.html

http://moltemplate.org/lammps_code/pair_lj_charmm_coul_charmm_inter.html
(Search for "soft-core".  Note: This feature requires you to recompile LAMMPS.)

### This model of HIV DNA lacks proteins or RNA

HIV is a retrovirus.
During the process of infection, the RNA genome is converted into DNA
in a process called reverse-transcription.  During this time, some of the
RNA will be degraded and will exit the virus.
This process is believed to occur inside the conical-shaped capsid shell.
(shown as grey spheres in the figure).
Some RNA may remain after the process is finished.
However RNA was not included in these simulations.

In the real HIV virus, the mature capsid is filled with nucleocapsid and
integrase proteins (among others) which are attracted to DNA and RNA.
These proteins will dramatically affect the conformation of the DNA and RNA
within the capsid, but they were also not included in this example.

I realize this is a complicated example.
Hopefully in spite of that, this example is useful.

Andrew, 2020-12-30
