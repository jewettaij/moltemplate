This directory contains examples of all-atom simulations using the DREIDING
force field.

As of 2020-10-19, this software has been tested for consistency between LAMMPS
and Materials Studio using only two simple molecules (butane and ethylene).

This software is experimental and the force-fields and equilbration protocols
have not been tested carefully.  There is no gaurantee that simulations
prepared using moltemplate will reproduce the behavior of Materials Studio
(or other MD codes).

Please report any bugs you find.

### *The atomic charges in these examples are not correct*

Some force-fields (like COMPASS, and moltemplate's version of OPLSAA) include
rules for assigning partial charges to atoms.  Most force fields, including
DREIDING do not.  So DREIDING users will have to obtain atomic charges by
some other means, perhaps by using 3rd-party tools.

WARNING:
*In these examples, I obtained partial charge estimates
from the OPLSAA parameter file located
[here](http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm).)*
***PLEASE DO NOT DO THIS!***
This is not how the DREIDING force field should be used.
It will probably not result in accurate behavior.

Alternatively LAMMPS'
[fix qeq/point](https://lammps.sandia.gov/doc/fix_qeq.html)
feature can be used to assign partial charges.
**QEq** is one of the built-in methods
[used by Gaussian to assign partial charges to atoms](https://gaussian.com/mm/)
in molecules that are modeled using the DREIDING force field.
If this fix is run infrequently (or if it is run only once at the beginning
of the simulation), then it should not slow the simulation down significantly.

In addition to the QEq method, the
[Gasteiger method](https://doi.org/10.1016/0040-4020(80)80168-2)
(suggested in the
[original DREIDING](https://doi.org/10.1021/j100389a010) paper)
[is also popular.](https://bioviacommunity.force.com/Communities_Topics?id=90638000000Gw1YAAS&name=Materials+Studio&tid=09a500000004QfAAAU#/feedtype=SINGLE_QUESTION_DETAIL)


### Improper angles

The style of improper interaction used by DREIDING force fields depends on an
angle which depends on the order of the atoms surrounding the central atom.
When multiple atoms have the same type, this creates ambiguity in atom order.
Since there is no guarantee that moltemplate will choose the same atom order
as other molecule builders this can lead to small unavoidable discrepancies
in energies and forces computed by LAMMPS and other simulation programs
(like Materials Studio).  But their effect should be neglegible.
(Please let us know if this is not the case.)
