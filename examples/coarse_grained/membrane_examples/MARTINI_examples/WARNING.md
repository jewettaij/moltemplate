WARNING regarding the MARTINI force field in moltemplate
====================

The MARTINI lipid examples have been implemented on two different
occasions by different people.

1) The current implementation was created by David Stelter and
Pieter J. in ’t Veld, and it uses v2.0 of the MARTINI force field.
(Note that prior to 2021-2-22, this version had incorrect electrostatics
and would produce incorrect results.)

2) An older implementation was created by Saeed M Bashusqeh in 2013.
(These files are in the "OLD_VERSION_wrong_number_of_angles" directory.
This version had too many angle interactions in the DPPC lipid.
Please avoid this version.  I will probably remove it eventually.)

## Disclaimer

I mention this in order to point out that previous iterations
of the MARTINI force field in moltemplate were incorrect.
There is no guarantee that the current version is correct.

The moltemplate developers have historically relied on input from
users to spot bugs in force fields.  There were several mistakes in our
implementations of the OPLSAA, GAFF2, COMPASS, and DREIDING force fields.
However these force fields were popular enough that users have spotted
these mistakes and reported them to us.  Unfortunately, coarse grained
force fields like MARTINI are less popular, and we haven't received
feedback yet for MARTINI.

While I would like to get MARTINI validated and working in moltemplate,
it is currently safer to use the official version of MARTINI
which is available [here](http://cgmartini.nl).

Andrew 2021-2-22
