3bp2p model
====================

This example demonstrates how to prepare a simulation of double-stranded (ds)
DNA using the coarse-grained "3bp2p" (3 base pairs = 2 particles) model.

### Images

![](https://raw.githubusercontent.com/jewettaij/dlpdb/master/examples/dna_example/statistics_keeping_every_3rd_base_pair/dsDNA_3to1_C3p.png)

*(additional images are shown below)*

During the simulation, the forces acting on the particles in the coarse-grained model are calculated from equations containing these degrees of freedom: r1, r2, θ1, θ2, φb, φt.  See the [dna_forcefield.lt](moltemplate_files/dna_forcefield.lt) file for details.  (The method used to generate these force field parameters is explained [here](moltemplate_files/deriving_force_field_parameters).)


### Features

1) This model has a double-helical conformation with major and minor grooves.
2) The polymer can be bent arbitrarily far without kinking or numerical instability.  This allows the use of comparatively large timesteps.
3) The masses were adjusted to make it possible to set the timestep to "1.0".  (*Note:* Use of a common timestep makes it easier to combine different coarse grained models together in the same simulation.  Using heavier masses can compensate for stiff bonds or bond-angles.  Using unrealistic masses is justified if the processes-of-interest occur on time-scales far larger than the vibrational spectra of individual bonds or angles in the coarse grained model.  At these time scales, the rates of diffusive motion can be controlled by adjusting the Langevin parameters and solvent properties, and are not effected by the mass of the particles of the model.)

Geometrical and mechanical properties of this DNA model (simulated at T=300K) are compared with experiment below:

| property                                       |  this model  | experimenal | ref |
|------------------------------------------------|--------------|-------------|-----|
| persistence length:                            |  48.7101 nm  |   46-48 nm  |  1  |
| torsional persistence length:                  | 109.9752 nm  | 100-110 nm  |  1  |
| helical twist angle between monomers (3bp):    | 102.7797 deg |  see below  |     |
| → helical twist angle between base-pairs:     |  34.2599 deg |   34.3 deg  |  2  |
| distance along the axis between monomers (3bp):|  0.98293 nm  |  see below  |     |
| → distance along the axis between base-pairs: |  0.32764 nm  |   0.332 nm  |  3  |

##    Prerequisites

LAMMPS must be compiled with the "MOLECULE" AND "EXTRA-MOLECULE"
packages enabled.  If LAMMPS generates the following error:
"dihedral_style spherical: Unknown dihedral style", then you must follow
[these instructions](https://docs.lammps.org/Build_package.html),
and recompile LAMMPS.

##    WARNING

These files (originally uploaded on 2019-12-10) contain many comments
(beginning with "#") which are probably misleading and no longer relevant.
I will try to clean up these files over time.
Please let me know if anything doesn't work.

##    Instructions

Instructions on how to build LAMMPS input files and 
run a short simulation are provided below

The following file contain instructions explaining how to generate
the curve that you want the DNA polymer to follow.
(You can also run it as an executable.)

   ./STEP_1_generate_initial_path.sh

The next file explains how to convert this curve into a moltemplate file, and
how to run moltemplate on that file. (You can also run it as an executable.)

   ./STEP_2_generate_LAMMPS_files.sh

Finally, to run the LAMMPS simulation follow the instructions in this file:
STEP_3_run_sim.sh
You will have to edit the file to specify the name of the LAMMPS binary
you are using (for example, "lmp_ubuntu"), and the number of processors

## References

1) Kriegel, F., Ermann, N., Forbes, R., Dulin, D., Dekker, N.H. and Lipfert, J., Nucleic Acids Research (2017), Vol. 45, No. 10, pp.5920-5929

2) Sinden, Richard R. (1994-01-15). "DNA structure and function" (1st ed.). Academic Press. p. 398. ISBN 0-12-645750-6

3) Theodore J. Zwang, Sylvia Hurlimann, Michael G. Hill, and Jacque-line K. Barton, J Am Chem Soc (2016), 138(48):15551-15554, "Helix-dependent Spin Filtering through the DNA Duplex"

## Examples

### Bacterial bottlebrush nucleoid model

![](http://moltemplate.org/images/misc/plectenomic_bottlebrush_off_axis_bbk.jpg)

This is a crude early model of a bacterial chromosome created by twisting a linear polymer under weak tension.  (The degree of supercoiling is unrealistically high in this example.)


### Confined supercoiled DNA simulation

![](http://moltemplate.org/images/prokaryotes/confined_supercoiling_2016-1-13.jpg)

This is a simulation of twisted DNA confined at concentrations mimicking the conditions inside *Bdellovibrio bacteriovorus* nucleoid.  *(Note: This is NOT a realistic model of a bacterial nucleoid.  This simulation has not been run long enough to reach equilibrium.  This model is 1/8th the size of the bacteria.  The degree of supercoiling is probably unrealistically high.)*/

### Bacterial bottlebrush nucleoid model (with crude NAPS)

![](http://moltemplate.org/images/prokaryotes/condensation_supercoiling+proteinglue_2016-3-26.jpg)

This is a crude early model of a small section of the *Caulobacter crescentus* nucleoid, simulated in the presence of small proteins (white spheres) which non-specifically stick to DNA.  This is not a realistic model.  (*Note: There's a video of simulation* [***here***](https://www.youtube.com/watch?v=A_ER8ztxl5I).)

### Wrapping DNA along a random space-filling curve

![DNA model](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/moltemplate_usage/CG_dsDNA_gold_turquoise.gif)
![](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/plus.png)
![space-filling curve](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/hamiltonian_paths_16x16x16.gif)
![](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/rightarrow.png)
![DNA wrapped along curve](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/moltemplate_usage/wrap_CG_dsDNA_around_a_curve_from_ndmansfield_LLR.png)

Simple random curves are often used as a starting point for simulations of DNA in confined conditions. (The process of DNA relaxation and unknotting via topoisomerases can be simulated later.)

*(Note: This simulation was prepared using [genpoly_lt.py](../../../../../../../doc/doc_genpoly_lt.md), [ndmansfield](https://github.com/jewettaij/ndmansfield), and [interpolate_curve.py](../../../../../../../doc//doc_interpolate_curve.md).)*
