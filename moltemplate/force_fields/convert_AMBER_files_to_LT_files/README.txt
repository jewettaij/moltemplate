## NOTE

The scripts needed to convert AMBER force-field files (eg "gaff.dat", *.frcmod)
into moltemplate (.LT) format have moved.

All of this code has been merged into a single file named "amber2lpt.py"
located in the "../../" directory (where most of the .py files are located).

Usage example 1:
```
amber2lt.py --in gaff.dat --out gaff.lt --name GAFF
```
Usage example 2:
```amber2lt.py --in DABPA_AC.frcmod --out forcefield_dabpa_ac.lt --name DABPA_AC_FF
```
