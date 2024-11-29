oplsaa2lt.py
==========================

**oplsaa2lt.py** converts a pair of BOSS input files (.par, .sb) into a moltemplate-compatible .LT file

# Typical Usage
If you have downloaded the latest BOSS parameter files
(eg "oplsaa.par", and "oplsaa.sb"),
you can convert them to a moltemplate file (eg. "oplsaa.lt") this way:
```
./oplsaa2lt.py --name OPLSAA --out oplsaa.lt \
  --par oplsaa.par \
  --sb  oplsaa.sb
```
In this example, the two BOSS files
- `oplsaa.par` contains atom definitions, dihedral and improper interaction paramters.
- `oplsaa.sb` contains bond and angle parameters.
*(Example .par and .sb files are located in the ../oplsaa2023_original_format/ directory.  Those files were published in 2023.)*

This creates a moltemplate file named "oplsaa.lt" containing a force-field object named "OPLSAA"
```
OPLSAA {
  :   (Atom types and force-field parameters go here...)
}
```
Later on, users can use this file to define molecules (eg. "Ethylene")
using this syntax:
```
import "oplsaa.lt"

Ethylene inherits OPLSAA {
  # atom-id mol-id atom-type charge   X       Y       Z     # comment
  write('Data Atoms') {
    $atom:c1  $mol @atom:143  0.00 -0.6695   0.00000  0.000  #143<->"alkene C (H2-C=)"
    $atom:h11 $mol @atom:144  0.00 -1.23422 -0.85446  0.000  #144<->"alkene H (H-C=)"
    :
  }
  # BondID     AtomID1  AtomID2
  write('Data Bond List') {
    $bond:c1h1 $atom:c1 $atom:h11
    :
  }
} # Ethylene
```
See the examples in the force_field_OPLSAA2023/ subdirectory directory for details.


The "oplsaa.lt" file is normally generated using the published parameters from [this paper](https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.3c06602).  Unless otherwise specified in the header, the `oplsaa2lt.py` program is run using these arguments:
```
./oplsaa2lt.py --name OPLSAA --out oplsaa.lt \
  --par ../oplsaa2023_original_format/Jorgensen_et_al-2023-The_Journal_of_Physical_Chemistry_B.sup-2.par \
  --sb  ../oplsaa2023_original_format/Jorgensen_et_al-2023-The_Journal_of_Physical_Chemistry_B.sup-3.sb
```
