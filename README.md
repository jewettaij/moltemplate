[![Build Status](https://travis-ci.org/jewettaij/moltemplate.svg?branch=master)](./.travis.yml)
[![GitHub](https://img.shields.io/github/license/jewettaij/moltemplate)](./LICENSE.md)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/moltemplate)](https://pypistats.org/packages/moltemplate)
[![PyPI - Version](https://img.shields.io/pypi/v/moltemplate)](https://pypi.org/project/moltemplate/)
[![Website](https://img.shields.io/website?down_color=orange&down_message=moltemplate.org%20offline&up_color=green&up_message=online&url=https%3A%2F%2Fmoltemplate.org)](http://moltemplate.org)
[![GitHub repo size](https://img.shields.io/github/repo-size/jewettaij/moltemplate)]()



Moltemplate
===========

##  Description

[Moltemplate](http://moltemplate.org)
is a *general* cross-platform text-based molecule builder for
[**LAMMPS**](https://lammps.sandia.gov) and *(the TCL version of)*
[**ESPResSo**](http://espressomd.org).
Moltemplate was intended for building custom coarse-grained molecular models,
but it can be used to prepare realistic all-atom simulations as well.
It currently supports the
[**OPLSAA**](./examples/all_atom/force_field_OPLSAA),
[**OPLSUA**](./examples/all_atom/force_field_OPLSUA_united_atom),
[**LOPLSAA**(2015)](./examples/all_atom/force_field_OPLSAA/hexadecane),
[**AMBER**(GAFF,GAFF2)](./examples/all_atom/force_field_AMBER),
[**DREIDING**](./examples/all_atom/force_field_DREIDING),
[**COMPASS**](./examples/all_atom/force_field_COMPASS),
[**TraPPE**(1998)](./examples/coarse_grained/solvent_models/manybodywaterMW+hydrocarbonsTraPPE)
force fields,
the
[**ATB**](https://atb.uq.edu.au) molecule database,
and the
[**MOLC**](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp04120f),
[**mW**](https://doi.org/10.1021/jp805227c),
[**ELBA**(water)](./examples/coarse_grained/solvent_models/ELBAwater%2Bmethanol),
[**oxDNA2**](https://dna.physics.ox.ac.uk/index.php/DNA_model_introduction),
and
[**EFF**](./examples/misc_examples/explicit_electrons/eff_CH4)
molecular models (and others).
Moltemplate is interoperable with
[**ATB**](https://atb.uq.edu.au),
[**VMD/topotools**](https://www.ks.uiuc.edu/Research/vmd),
[**PACKMOL**](http://m3g.iqm.unicamp.br/packmol/home.shtml),
[**Open Babel**](https://open-babel.readthedocs.io/en/latest/FileFormats/The_LAMMPS_data_format.html),
[**EMC**](http://montecarlo.sourceforge.net/),
[**CellPACK**](http://www.cellpack.org),
[**LigParGen**](http://moltemplate.org/doc/moltemplate_talk_2019-8-15.pdf#page=190),
[**Vipster**](https://sgsaenger.github.io/vipster),
[**struc2lammpsdf**](https://nanohub.org/resources/struc2lammpsdf),
and any other program that reads or generates LAMMPS data (.lmpdat) files.
(New force fields and examples are added continually by users.)


This repository contains 3 folders:
-[moltemplate](./moltemplate/): source code and force fields
-[doc](./doc/): documentation for moltemplate, ltemplify, genpoly_lt, etc..
-[examples](./examples/): examples built with moltemplate

## Installation Instructions [Here](INSTALL.md)

### Documentation

The best way to learn how to use moltemplate is to find an example
which is similar to the system that you wish to simulate and modify it.
This repository includes approximately 50 [examples](./examples).
Some of the moltemplate examples are also demonstrated (with pictures)
[here](http://moltemplate.org/visual_examples.html).

All moltemplate users should probably read chapter 4 of the
[reference manual](./doc/moltemplate_manual.pdf)
*(It's only a few pages long.  The first 3 chapters are optional.)*
In addition, there are also several
[talks/tutorials](http://moltemplate.org/doc/talks.html)
online.


### Limitations for preparing all-atom simulations

Moltemplate [does *not* choose atom types automatically ("atom typing")](http://moltemplate.org/force_field_recommendations.html),
and currently cannot be used to build all-atom proteins from scratch.
3rd-party tools may be needed to calculate atomic partial charges accurately.
(For non-trivial molecules, users are encouraged to use the
[ATB](https://atb.uq.edu.au) database,
or use the DREIDING force field which has simple [atom type rules](./doc/DREIDING_Label_Manual.pdf),
or use a 3rd-party molecule-builder
and convert the resulting files to LAMMPS format using
[OpenBabel](https://open-babel.readthedocs.io/en/latest/FileFormats/The_LAMMPS_data_format.html "Convert 3rd party sim files to LAMMPS DATA format")
followed by [ltemplify.py](./doc/doc_ltemplify.md "Convert LAMMPS DATA to Moltemplate format"), and then use moltemplate.)



## Typical usage
```
    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt
```

## Web page

Additional suggestions and supporting code can be found at:

http://www.moltemplate.org


## Requirements

Moltemplate requires the Bourne-shell, and a recent version of python
(2.7, 3.0 or higher), and can run on OS X, linux, or windows.
(...if a suitable shell environment has been installed.  See below.)

The *numpy* python module is also required.

## License

With the exception of one file
([ttree_lex.py](./moltemplate/ttree_lex.py)),
moltemplate is available under the terms of the [MIT license](LICENSE.md).

The remaining file, ([ttree_lex.py](./moltemplate/ttree_lex.py)),
is a modified version of the 
[shlex.py](https://docs.python.org/3/library/shlex.html) module,
which was released using the
[PSF license](https://docs.python.org/3/license.html).
Hence [ttree_lex.py](./moltemplate/ttree_lex.py) must also use this license.
(*The PSF is not a copyleft license.
It is similar to the BSD and MIT licenses and
[is compatible with the the GPL license](https://docs.python.org/3/license.html).)*


## Citation

If you find this program useful, please cite:

*"Moltemplate: A Tool for Coarse-Grained Modeling of Complex Biological Matter and Soft Condensed Matter Physics", J.Mol.Biol., (2021), **in press**, Jewett AI, Stelter D, Lambert J, Saladi SM, Roscioni OM; Ricci M, Autin L, Maritan M, Bashusqeh SM, Keyes T, Dame RT; Shea J-E, Jensen GJ, Goodsell DS*
