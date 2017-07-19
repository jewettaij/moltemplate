This is an initial test of the "molecular cellular automata" extensions for
LAMMPS and moltemplate.

Motivation:
Frequently in biology proteins undergo changes when interacting with other
proteins.  Chains of such transitions can lead to complex behavoir, such
as the "dynamic instability" found in microtubules and the ParM system.

I am using custom versions of "fix bond/break" and "fix bond/create"
(currently named "fix bond/change" and "fix bond/new") which make this possible.

Youtube video of this system:

https://www.youtube.com/watch?v=EEbt07vZHew
(explanation and color legend below)

The system simulated here was inspired by the ParM
bacterial chromosome segregation system. (ParM's assembly/dissasembly
mechanism is also very similar to the mechanism of microtubule growth/collapse
in eukarya.) In this coarse grained simulation, each particle represents a
nucleotide-bound-ParM protein whose color indicates whether it is bound to ATP
(light-green,green) or ADP (light-blue,blue). (The ATP molecules are not
modeled explicitly.) Initially, free ATP-bound-ParM aggregate into short
polymers and grow by accumulating at the end (green). Over time, the
ATP-bound-ParM monomers hydrolyze, converting to ADP-bound-ParM monomers (blue).
(ATP-bound-ParM proteins do not hydrolyze ATP quickly unless part of a polymer.)
ADP-bound-ParM is unstable and will rapidly unbind from the end of a polymer,
causing the sudden depolymerization events seen at 0:11, 0:21, and 0:41.
Later these free ADP-bound-ParM proteins revert back to ATP-bound-ParM,
and the cycle repeats itself. This simulation can be built with
approximately ~100 lines of moltemplate & LAMMPS code and it captures the
cooperative process of polymer (trimer) nucleation, growth, collapse, and ATP
recycling.
Note: This is a tech demonstration. The parameters here have not yet
been tuned to reproduce realistic kinetics of ParM filament assembly in
bacteria. (Real ParM filaments are helical filaments formed by ParM dimers,
not monomers, and their growth is not bidirectional as shown here, but is
stimulated in-vivo at one end by ParRC.)

legend:
TF   <-->   ParM+ATP free  (light green in the youtube video)
TB   <-->   ParM+ATP bound  (green in the youtube video)
TE   <-->   ParM+ATP endcap  (dark green in the youtube video)
DF   <-->   ParM+ADP free   (light blue in the youtube video)
DB   <-->   ParM+ADP bound  (blue in the youtube video)
DE   <-->   ParM+ADP endcap  (dark blue in the youtube video)
