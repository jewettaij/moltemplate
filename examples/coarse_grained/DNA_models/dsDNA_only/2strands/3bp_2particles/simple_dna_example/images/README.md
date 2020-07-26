images
========

To keep the repository size small, pictures of this coarse-grained model were
omitted.  However some images of this DNA model can be found online.
(Links below.)


### Model Details:

![](https://raw.githubusercontent.com/jewettaij/dlpdb/master/examples/dna_example/statistics_keeping_every_3rd_base_pair/dsDNA_3to1_C3p.png)

(Note: During the simulation, the forces acting on the particles in the 
coarse-grained model are calculated from equations containing these
degrees of freedom: r1, r2, θ1, θ2, φb, φt.
Unfortunatley, the names and definitions of the θ and φ angles may have changed.
*I really should update this figure.*
***Ghastly*** details can be found [here](../moltemplate_files).)


### Bacterial bottlebrush nucleoid model

![](http://moltemplate.org/images/misc/plectenomic_bottlebrush_off_axis_bbk.jpg)

This is a crude early model of a bacterial chromosome created by twisting a linear polymer under weak tension.  (The degree of supercoiling is unrealistically high in this example.)

### Confined supercoiled DNA simulation

![](http://moltemplate.org/images/prokaryotes/confined_supercoiling_2016-1-13.jpg)

This is a simulation of twisted DNA confined at concentrations mimicking the
conditions inside *Bdellovibrio bacteriovorus* nucleoid.
*(Note: This is NOT a realistic model of a bacterial nucleoid.
 This simulation has not been run long enough to reach equilibrium.
 This model is 1/8th the size of the bacteria.
 The degree of supercoiling is probably unrealistically high.)*/


### Bacterial bottlebrush nucleoid model (with crude NAPS)

![](http://moltemplate.org/images/prokaryotes/condensation_supercoiling+proteinglue_2016-3-26.jpg)

This is a crude early model of a small section of the
*Caulobacter crescentus* nucleoid, simulated in the presence of small
proteins (white spheres) which non-specifically stick to DNA.
This is not a realistic model.
(*Note: There's a video of simulation* [***here***](https://www.youtube.com/watch?v=A_ER8ztxl5I).)

### Wrapping DNA along a random space-filling curve

![DNA model](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/moltemplate_usage/CG_dsDNA_gold_turquoise.gif)
![](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/plus.png)
![space-filling curve](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/hamiltonian_paths_16x16x16.gif)
![](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/rightarrow.png)
![DNA wrapped along curve](https://raw.githubusercontent.com/jewettaij/ndmansfield/master/doc/images/moltemplate_usage/wrap_CG_dsDNA_around_a_curve_from_ndmansfield_LLR.png)

Simple random curves are often used as a starting point for simulations
of DNA in confined conditions.
(The process of DNA relaxation and unknotting via topoisomerases
 can be simulated later.)

*(Note: This simulation was prepared using [genpoly_lt.py](https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_lt.md), [ndmansfield](https://github.com/jewettaij/ndmansfield), and [interpolate_curve.py](https://github.com/jewettaij/moltemplate/blob/master/doc/doc_interpolate_curve.md).)*

