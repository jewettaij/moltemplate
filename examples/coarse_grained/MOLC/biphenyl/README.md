## WARNING

This example uses several new LAMMPS force field styles which were
implemented by the authors of the MOLC paper (see below), and are
not yet available in LAMMPS (as of 2020-6-01).
Hopefully by the time you are reading this, that has changed.
The syntax of the commands used to specify force field parameters
may change as the MOLC code evolves.  I will try to keep this
example up to date.  Please report any problems to the github
issue tracker  https://github.com/jewettaij/moltemplate/issues   
or to my email (jewett dot aij  at  gmail dot com)

## Description

MOLC is a new strategy for generating coarse grained versions of
all-atom molecules using bonded ellipsoids. In this example, we have
used the MOLC strategy to study a sample of 1000 biphenyl molecules.
First, an all-atom representation of the molecule (biphenyl) is created,
for example using the ATB molecule server
   https://atb.uq.edu.au/
Then it is converted to a MOLC coarse grained model using the procedure
outlined in this paper:

*"MOLC. A reversible coarse grained approach using anisotropic beads for
the modelling of organic functional materials"
Ricci M, Roscioni OM, Querciagrossa L, Zannoni C, Phys.Chem.Chem.Phys.
(2019), 21(47):26195-26211*

The MOLC version of the biphenyl molecule is represented by only two particles.
The intra-molecular interactions are represented by a single bond.
The inter-molecular interactions are described with a short-range Gay-Berne
potential and a long-range Coulomb potential. The sample was condensed from
a cubic grid of molecules, using 3D periodic boundary conditions and a time
step of 10 fs. The temperature was controlled with a Langevin thermostat with
coupling constant (damping time) of 1 ps. A trajectory of 40 ns was produced
in the canonical ensemble at T = 300 K, P = 1 atm and a timestep of 20 fs.
A coupling constant of 10 ps was used for the Nosé-Hoover barostat.
The density of the resulting sample is 1.040(4) g/cm3. The density
and radial distribution function shows an excellent agreement with
an all-atom reference model.
The distribution of the dihedral angles between the phenyl rings is
consistent between the coarse-grained MOLC model and the all-atom model,
with the atomistic model having a tilt angle of 51 degrees and the MOLC
model an angle of 66 degrees. 
