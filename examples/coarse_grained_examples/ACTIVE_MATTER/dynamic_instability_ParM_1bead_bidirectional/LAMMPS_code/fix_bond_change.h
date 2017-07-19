/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#define ALLOW_BOND_MODIFY_SINGLE_ATOMS
#define ALLOW_BOND_MODIFY_BREAK_BONDS

#ifdef FIX_CLASS

FixStyle(bond/change,FixBondChange)

#else

#ifndef LMP_FIX_BOND_CHANGE_H
#define LMP_FIX_BOND_CHANGE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondChange : public Fix {
 public:
  FixBondChange(class LAMMPS *, int, char **);
  ~FixBondChange();
  int setmask();
  void init();
  void post_integrate();
  void post_integrate_respa(int,int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me,nprocs;
  //int btype;
  int seed;
  int angleflag,dihedralflag,improperflag;
  bigint nevery_delay;
  bigint lastcheck;

  // Variables that specify the condition(s) for a transition to occur
  // as well as the final state after the transition
  // If multiple conditional requirements are requested, all must be satisfied.
  double bond_rminsq, bond_rmaxsq;        // select bonds by distance
  int atype1old_lo, atype1old_hi, atype2old_lo, atype2old_hi;//select atom types
  int atype1new, atype2new;                                  //new atom types
  double qold1_lo, qold1_hi, qold2_lo, qold2_hi; // select atoms by charge
  double qnew1, qnew2;                           // new charges after transition
  int btypeold_lo, btypeold_hi;           // select bond type before transition
  int btypenew;                           // new bond type after transition
  double prob_transition;
  int conserve_ke_flag;           // 1 = conserve ke, 0 = do not conserve ke
  bool check_bonded_atoms;
  bool change_bond_properties;
  bool prioritize_long_bonds;
  bool break_bond;
  int ignore_tags_flag; //0 = restriction for atoms i,j to satisfy tag[i]<tag[j]
                        //1 = no restrictions based on atom tag
  //bool iterate_until_stationary;  (not implemented yet)

  #ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS
  int atypeold_lo, atypeold_hi, atypenew; // change atom type regardless of bond
  double qold_lo, qold_hi, qnew;          // change charge regardless of bond
  bool check_individual_atoms;
  #endif


  int breakcount,breakcounttotal;
  int nmax;
  tagint *partner,*finalpartner;
  double *distsq,*probability;

  int nbreak,maxbreak;
  tagint **broken;

  tagint *copy;

  class RanMars *random;
  int nlevels_respa;

  int commflag;
  int nbroken;
  int nangles,ndihedrals,nimpropers;

  bool SatisfiesAtomProperties(int, int);
  void check_ghosts();
  void update_topology();
  void break_angles(int, tagint, tagint);
  void break_dihedrals(int, tagint, tagint);
  void break_impropers(int, tagint, tagint);
  void rebuild_special_one(int);
  int dedup(int, int, tagint *);

  // DEBUG

  void print_bb();
  void print_copy(const char *, tagint, int, int, int, int *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid bond type in fix bond/change command

Self-explanatory.

E: Cannot use fix bond/change with non-molecular systems

Only systems with bonds that can be changed can be used.  Atom_style
template does not qualify.

E: Cannot yet use fix bond/change with this improper style

This is a current restriction in LAMMPS.

E: Fix bond/change needs ghost atoms from further away

This is because the fix needs to walk bonds to a certain distance to
acquire needed info, The comm_modify cutoff command can be used to
extend the communication range.

*/
