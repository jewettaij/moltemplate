Note: By default, the system.data and system.in.settings files contain
extra information for atoms defined in DREIDING ("dreiding.lt") which you
are not using in this simulation.
This is harmless, but if you to delete this information from your
system.in.settings and system.in.data files:


1) run this script:

cleanup_moltemplate.sh

(Note: Removing unecessary atom types will make it easier to visualize the
       simulation in VMD.)


The user will need to modify the system.in.init and system.in.settings files.
This is because DREIDING uses a special hydrogen bonding term however, furan
has no hydrogen bonds.  (Consequently running "cleanup_moltemplate.sh" will
remove the definitions of atoms which can make hydrogen bonds, and this will
cause LAMMPS to halt with an error message.)


2) Edit the "system.in.init" file and replace:

pair_style hybrid/overlay hbond/dreiding/lj 4 6 6.5 90 lj/cut/coul/long 10.0

with

pair_style lj/cut/coul/long 10.0

If all of the atoms in your molecule contain no charge (which is unlikely), then
use this line instead:

pair_style lj/cut 10.0

...and delete this line:

kspace_style ewald 0.0001

3) Edit the "system.in.settings" file, find the lines that begin with
"pair_coeff", and delete the phrase "lj/cut/coul/long" from each of these lines.

