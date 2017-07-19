Note:

fix_bond_new.cpp is my customized version of fix_bond_create.cpp

fix_bond_change.cpp is my customized version of fix_bond_change.cpp
(Most of the work so far has been getting fix_bond_change.cpp to work.)

fix bond/change allows atom types to change using this syntax:

  atoms i_old j_old -> i_new j_new 

(I intend to change the syntax of the fix_bond_new.cpp 
command so that it matches the syntax used by fix_bond_change.cpp,
but I haven't got around to that yet.)

Breaking bonds is optional in fox bond/change.  To force a bond to break, use:

  bond BREAK

Complex changes in the state of the system can be caused by stringing
a chain of fix bond/new and fix bond/change commands together in succession.

The optional "ndelay" argument allows you to delay each fix by "ndelay"
timesteps to make sure that a chain of fixes are evaluated in the correct order.
Typically "ndelay" is less than the interval between fix attempts ("nevery").

I also intend to introduce a new fix fix_atom_change.cpp which
can stochastically change single atoms one at a time (not in pairs).
(Right now this functionality is handled inside "fix_bond_change.cpp".)

