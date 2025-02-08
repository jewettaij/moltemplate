[![CircleCI](https://circleci.com/gh/jewettaij/moltemplate.svg?style=svg)](https://circleci.com/gh/jewettaij/moltemplate)
[![CodeQL](https://github.com/jewettaij/moltemplate/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/jewettaij/moltemplate/actions/workflows/codeql-analysis.yml)
[![Requires Python 3.4+](https://img.shields.io/badge/Python-3.4+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI - Version](https://img.shields.io/pypi/v/moltemplate)](https://pypi.org/project/moltemplate/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/moltemplate?logo=pypi&logoColor=white&color=blue&label=PyPI)](https://pypistats.org/packages/moltemplate)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/moltemplate?logo=condaforge&color=blue&label=Conda)](https://anaconda.org/conda-forge/moltemplate)
[![Paper](https://img.shields.io/badge/J.MolBiol-2021.5.28-blue?logo=elsevier&logoColor=white)](https://doi.org/10.1016/j.jmb.2021.166841)
[![Website](https://img.shields.io/website?down_color=orange&down_message=moltemplate.org%20offline&up_color=green&up_message=online&url=https%3A%2F%2Fmoltemplate.org)](https://moltemplate.org)





Moltemplate
===========

##  Description

[Moltemplate](https://moltemplate.org)
is a *general* cross-platform text-based molecule builder for
[**LAMMPS**](https://lammps.sandia.gov).
Moltemplate was intended for building custom coarse-grained molecular models,
but it can be used to prepare realistic all-atom simulations as well.
It currently supports the
[**ATB**](https://atb.uq.edu.au) molecule database, the
[**OPLSAA**(2023)](https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.3c06602),
[**OPLSUA**(2008)](./examples/all_atom/legacy_force_field_examples/force_field_OPLSUA_united_atom_2008),
[**LOPLS**(2015)](https://doi.org/10.1021/acs.jpcb.5b08569),
[**AMBER**(GAFF,GAFF2)](./examples/all_atom/force_field_AMBER),
[**DREIDING**](./examples/all_atom/force_field_DREIDING),
[**COMPASS**](./examples/all_atom/force_field_COMPASS),
[**TraPPE**(1998)](./examples/coarse_grained/solvent_models/manybodywaterMW+hydrocarbonsTraPPE),
[**EFF**](./examples/misc_examples/explicit_electrons/eff_CH4),
force fields,
and the
[**MOLC**](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp04120f),
[**mW**](https://doi.org/10.1021/jp805227c),
[**ELBA**(water)](./examples/coarse_grained/solvent_models/ELBAwater%2Bmethanol),
and
[**oxDNA2**](https://dna.physics.ox.ac.uk/index.php/DNA_model_introduction),
molecular models (and others).
(New force fields and examples are added continually by users.)
Moltemplate is inter-operable with
[**ATB**](https://atb.uq.edu.au),
[**Open Babel**](https://open-babel.readthedocs.io/en/latest/FileFormats/The_LAMMPS_data_format.html),
[**AmberTools**](https://ambermd.org/AmberTools.php),
[**LigParGen**](https://moltemplate.org/doc/moltemplate_talk_2019-8-15.pdf#page=191),
[**RED-server**](https://upjv.q4md-forcefieldtools.org),
[**VMD**](https://www.ks.uiuc.edu/Research/vmd),
[**topotools**](https://sites.google.com/site/akohlmey/software/topotools/tutorial-introduction),
[**PACKMOL**](http://m3g.iqm.unicamp.br/packmol/home.shtml),
[**EMC**](http://montecarlo.sourceforge.net/),
[**CellPACK**](http://www.cellpack.org),
[**Vipster**](https://sgsaenger.github.io/vipster),
[**struc2lammpsdf**](https://nanohub.org/resources/struc2lammpsdf),
and any other program that generates [**MOL2**](https://github.com/UnixJunkie/mol2-file-format-spec/blob/master/mol2.pdf) or [**LAMMPS DATA**](https://docs.lammps.org/2001/data_format.html) (.lmpdat) files *(by using the [mol22lt.py](https://github.com/jewettaij/moltemplate/blob/master/doc/doc_mol22lt.md) and [ltemplify.py](https://github.com/jewettaij/moltemplate/blob/master/doc/doc_ltemplify.md) file converters)*.


This repository contains 3 folders:
- [moltemplate](./moltemplate/): source code and force fields
- [doc](./doc/): documentation for moltemplate.sh, ltemplify.py, genpoly_lt.py, etc..
- [examples](./examples/): examples built with moltemplate


### Documentation

The best way to learn how to use moltemplate is to find an example
which is similar to the system that you wish to simulate and modify it.
This repository includes approximately 50 [examples](./examples).
Some of the moltemplate examples are also demonstrated (with pictures)
[here](https://moltemplate.org/visual_examples.html).

All moltemplate users should read chapter 4 of the
[reference manual](https://moltemplate.org/doc/moltemplate_manual.pdf#section.4)
*(It's only 6 pages long. The first 3 chapters can be skipped.)*
Chapters 6, 7, and 5.1-5.3 also also recommended.
*(13 pages.)*

In addition, there are also several
[talks/tutorials](https://moltemplate.org/doc/talks.html)
online.


### Strategies for preparing all-atom simulations

Moltemplate was originally designed for preparing coarse-grained simulations.
Fully atomistic simulations usually require choosing
force-field specific atom types ("atom typing"),
and the calculation of atomic charge.
Moltemplate
[does not do this](https://moltemplate.org/force_field_recommendations.html).
Moltemplate also does not repair incomplete PDB files.
To get around these limitations:

* Users can use use *3rd-party molecule-builder tools* such as:
the [ATB](https://atb.uq.edu.au) database (LT format),
[AmberTools](https://ambermd.org/AmberTools.php) (MOL2 format),
[LigParGen](https://moltemplate.org/doc/moltemplate_talk_2019-8-15.pdf#page=190) (DATA format),
[OpenBabel](https://open-babel.readthedocs.io/en/latest/FileFormats/The_LAMMPS_data_format.html "Convert 3rd party sim files to LAMMPS DATA format"), and
[EMC](http://montecarlo.sourceforge.net/) (DATA format).
These tools can generate MOL2 or LAMMPS-DATA files
which can be converted into moltemplate format (LT format) using
[mol22lt.py](./doc/doc_mol22lt.md#mol22ltpy) or
[ltemplify.py](./doc/doc_ltemplify.md#ltemplifypy).
*(If you want to try and select the atom types
manually instead of using 3rd-party tools, read
[these suggestions](https://moltemplate.org/force_field_recommendations.html).)*
Alternatively, users can use the
[DREIDING force field](./examples/all_atom/force_field_DREIDING)
which has relatively simple
[atom type rules](./doc/DREIDING_Label_Manual.pdf).
* For atomic partial charges
the [RED-server](https://upjv.q4md-forcefieldtools.org)
is useful (which can generate a MOL2 file),
*in addition to the tools listed above*.
Alternatively, users can try using LAMMPS' built-in
[fix qeq](https://docs.lammps.org/fix_qeq.html) feature.



## Typical usage

```
    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt
```


## Installation Instructions

Moltemplate depends on other software to work
(such as BASH, python3, pip3, *or* conda).  Once these
dependencies have been met, installing moltemplate is relatively easy.
However many users find it difficult to install these prerequisites correctly.
A detailed installation guide is located [here](./INSTALL.md).


### Installation for conda users
If you have installed conda, you can download moltemplate using:
```
conda install conda-forge::moltemplate
```
This will omit the [moltemplate examples](./tree/master/examples)
and [documentation](./tree/master/doc).
You can browse examples [online](./tree/master/examples)
and/or download them using:
```
git clone https://github.com/jewettaij/moltemplate DESTINATION_DIRECTORY
```

### Installation for pip users

You can download moltemplate using:
```
git clone https://github.com/jewettaij/moltemplate DESTINATION_DIRECTORY
```
*(See below if you don't have git installed.)*
Then enter the directory where this README file is located:
```
cd DESTINATION_DIRECTORY
```
...and run the following command:
```
pip3 install .              # (or "pip", if that fails)
```
*(Note: In some environments, "pip3" is called "pip" instead.)*

If the command above fails (with both "pip" and "pip3"), then try this instead:
```
pip3 install . --user       # (or "pip", if that fails)
```
This will install moltemplate for a single user.
If you are on a shared computer and you want to install moltemplate
system-wide, then use:
```
sudo pip3 install .         # (or "pip", if that fails)
```
Later, you can uninstall moltemplate using:
```
pip3 uninstall moltemplate
# (use "pip" and/or prepend "sudo" if you did that earlier)
```
If this fails then read the
[installation troubleshooting guide](./INSTALL.md).
This guide will offer several different installation methods
and explain how to install the prerequisites needed.

***Note:***
*Alternatively, you can download and install moltemplate
using a single command:*
```
pip3 install git+https://github.com/jewettaij/moltemplate.git   # (or "pip", if that fails)
```
*...however this will omit all of the examples and documentation.*


## Web page

Additional suggestions and supporting code can be found at:

http://www.moltemplate.org


## Requirements

Moltemplate requires BASH and Python (3.4+), NumPy,
and can run on MacOS, linux, or windows.

*(Note: On MacOS, it may eventually become necessary to use 3rd-party tools
like "brew" to install BASH if apple removes BASH support in future updates.
Python 2.7 may also work, but you must edit your
[PATH](./INSTALL.md#Installation-method-2-Editing-bashrc)
to install moltemplate instead of using pip/pip3.
Pip no longer works with old versions of python.)*

To use LAMMPS and moltemplate, you will also need to install (and learn how to
use) a (unix-style) text editor.  (Word, Wordpad, and Notepad will not work.)
Popular graphical text editors
include **Atom**, **Sublime**, **Notepad++**, and **VSCode**.
Older, non-graphical programs include **vim**, **emacs**,
**nano**, **ne**, and **jove**.
(Apple's TextEdit can be used if you save the file as *plain text*.)


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

*"Moltemplate: A Tool for Coarse-Grained Modeling of Complex Biological Matter and Soft Condensed Matter Physics", J. Mol. Biol., 2021, 433(11):166841, Jewett AI, Stelter D, Lambert J, Saladi SM, Roscioni OM, Ricci M, Autin L, Maritan M, Bashusqeh SM, Keyes T, Dame RT, Shea J-E, Jensen GJ, Goodsell DS*
[https://doi.org/10.1016/j.jmb.2021.166841](https://doi.org/10.1016/j.jmb.2021.166841)
