#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
# License: MIT License  (See LICENSE.md)
# Copyright (c) 2020, Scripps Research
# All rights reserved.

man_page_text = """
Usage (example):

postprocess_coeffs.py ttree_assignments.txt < file.template

Where "file.template" contains contents similar to:

if atoms @atom:A @atom:B* @atom:C and 
   bond_type[1,2] == @bond:AB and
   bond_type[2,3] == @bond:BC and
   angle_type[1,2,3] == @angle:ABC and
   distance[1,3] < 3.0 and prob(0.5)
then
   atom_type[2] = @atom:B and
   atom_type[3] = @atom:D and
   delete_bond(1,2) and
   bond_type[1,3] = @bond:AD

if atom_type[1] == @atom:A and 
   atom_type[2] == @atom:B* and
   atom_type[3] == @atom:C and
   bond_type[1,2] == @bond:AB and
   bond_type[2,3] == @bond:BC and
   angle_type[1,2,3] == @angle:ABC and
   distance[1,3] < 3.0 and prob(0.5)
then
   atom_type[2] = @atom:B and
   atom_type[3] = @atom:D and
   delete_bond(1,2) and
   bond_type[1,3] = @bond:AD

if atom_type[1] == @atom:A and 
   atom_type[2] == @atom:B* and
   atom_type[3] == @atom:C and
   bond_type[1,2] == @bond:AB and
   bond_type[2,3] == @bond:BC and
   rmsd((1,2,3), ((1.3,5.2,1.2),(2.4,4.5,6.6),(0.01,1.5,9.55)) <= 3.2
   then
   delete_bond(2,3)
   delete_atom(3)

if atom_type[1] == @atom:A and 
   atom_type[2] == @atom:B* and
   atom_type[3] == @atom:C and
   bond_type[1,2] == @bond:AB and
   bond_type[2,3] == @bond:BC and
   rmsd((1,2,3), ((1.3,5.2,1.2),(2.4,4.5,6.6),(0.01,1.5,9.55)) <= 3.2
then
   coords((1.3,5.2,1.2),(2.4,4.5,6.6),(0.01,1.5,9.55),(-1.2,0.1,12))  #add atom
   and atom_type[4] = @atom:D
   and bond_type[3,4] = @bond:CD

"""

#OLD SYNTAX:
#
#if atoms((1,@atom:A), (2,@atom:B*), (3,@atom:C)) and &
#   bonds(((1,2),@bond:AB), ((2,3),@bond:BC)) and &
#   angles(((1,2,3),@angle:ABC)) and &
#   distance(1,3) < 3.0 and &
#   prob(0.5)
#then &
#   atoms( ,@atom:B,@atom:D) and &
#   bonds(((1,2),@bond:*), ((1,3),@bond:AD))
#
#if atoms(@atom:A, ,@atom:C) and &
#   bonds((1,2), (2,3)) and &
#   angles(((1,2,3),@angle:ABC)) and &
#   distance(1,3) < 3.0 and &
#   prob(0.5) and &
#then &
#   atoms(@atom:A, @atom:*,@atom:D) and &
#   bonds(((1,2),@bond:*), ((1,3),@bond:AD))






#Strategy
#1) Instantiate the lexer
#2) Use lex.ReadTemplate()   (and also StaticObj.CleanupReadTemplate() ?? no.)
#3) Augment the functionality of ttree_lex.SplitTemplateMulti() to optionally
#   enable it to not discard the delimiter used for splitting.  This way you
#   can see which delimiter appears where.  (This is optional and not default)
#4) Use SplitTemplateMulti(['if','then'])
#5) Use SplitTemplate('and')
#5) Use SplitTemplateMulti(['[', ']', '(', ')', ',', '==', '=', 'atom_type', 'bond_type', 'angle_type', 'dihedral_type', 'improper_type', 'distance', 'prob', 'rmsd', 'coords', 'delete_atom', 'delete_bond', 'delete_angle', 'delete_dihedral', 'dihedral_improper'])
#
# ------ the following features not needed for version 1.0: ------
#
#6) This next step is only needed for rmsd((),())  and coords():
#   Create a variant of SplitQuotedString() that splits according to both
#   parenthesis and commas and works with nested expressions.
#   Call it "SplitNestedQuotes()".  It accepts these arguments with these
#   default values:
#      delim = ','
#      open_paren = '(',
#      close_paren = ')',
#   It will split template lists of this form
#    ['bond((1,2),', VarRef('@/bond:AB'), '), ((2,3),', VarRef('@/bond:BC'),'))']
#        ... which is what ReadTemplate() will return when invoked on
#            'bond(((1,2),@/bond:AB), ((2,3),@/bond:BC))'
#   into something like this:
#    ['bond',
#     '(',
#     ['(1,2),', VarRef('@/bond:AB')],
#     ['(2,3),', VarRef('@/bond:BC')],
#     ')']
#   Note: This function only processes the outermost paren expression.
#         The '(1,2),' was not processed further.  Had it been, it would have
#         returned [['(', ['1','2'] ,')'], '']
#   
#7) Use SplitNestedQuotes() to find the arguments following the
#   rmsd() and coords() keywords.

import sys
import argparse
from collections import defaultdict
import re
import gc

try:
    from .ttree import ExtractFormattingCommands
    from .ttree_lex import *

except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree import ExtractFormattingCommands
    from ttree_lex import *


g_filename = __file__.split('/')[-1]
g_module_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str = '2020-10-24'
g_version_str = '0.0.1'
g_program_name = g_filename
#sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')


def ExtractVarName(text):
    """ Read a string like 'atom:A  '  or  '{/atom:A B/C/../D }ABC '
        and return ('','atom:A','  ')  or  ('{','atom:A B/C/../D ','}ABC')
        These are 3-tuples containing the portion of the text containing 
        only the variable's name (assumed to be within the text),
        ...in addition to the text on either side of the variable name.
    """
    i_begin = 0
    escape = '\''
    lparen = '{'
    rparen = '}'
    escaped = False
    commenters = '#'
    whitespace = ' \t\r\f\n'
    terminators = whitespace + commenters
    # Ideally, perhaps I should lookup these values from ttree_lex.TtreeLex to 
    # make sure I am being consistent, instead of re-defining them in this file.
    #while ((i_begin < len(text)) and
    #       (text[i_begin] in whitespace)):
    #    i_begin += 1
    in_paren = text[i_begin:i_begin+1] == lparen
    if in_paren:
        terminators = rparen
        i_begin += 1
    i_end = i_begin
    while ((i_end < len(text)) and
           (text[i_end] not in terminators)):
        i_end += 1
    return (text[0: i_begin],
            text[i_begin:i_end],
            text[i_end:])



def main():
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument('bindings_filename',
                        help='ttree assignments file name (usually "ttree_assignments.txt")')
        ap.add_argument('-t', '--template', required=False,
                        help='template text file (typically generated by moltemplate, and ending in ".template")')
        args = ap.parse_args()
        print(args.bindings_filename, args.template)
        exit(0)
        bindings_filename = ap.bindings_filename
        f = open(bindings_filename)
        atom_types = set([])
        bond_types = set([])
        angle_types = set([])
        dihedral_types = set([])
        improper_types = set([])

        #BasicUIReadBindingsStream(assignments, f, bindings_filename)

        # The line above is robust but it uses far too much memory.
        # This for loop below works for most cases.
        for line in f:
            #tokens = lines.strip().split()
            # like split but handles quotes
            tokens = SplitQuotedString(line.strip())
            if len(tokens) < 2:
                continue
            if tokens[0].find('@') != 0:
                continue
            if tokens[0][2:].find('atom') == 0:
                atom_types.add(tokens[0][1:])
            elif tokens[0][2:].find('bond') == 0:
                bond_types.add(tokens[0][1:])
            elif tokens[0][2:].find('angle') == 0:
                angle_types.add(tokens[0][1:])
            elif tokens[0][2:].find('dihedral') == 0:
                dihedral_types.add(tokens[0][1:])
            elif tokens[0][2:].find('improper') == 0:
                improper_types.add(tokens[0][1:])

        f.close()
        gc.collect()

        if ap.template is not None:
            templ_file = fopen(ap.template, 'r')
            lex = LineLex(ap.template)
        else:
            templ_file = sys.stdin
            lex = LineLex('__standard_input_for_postprocess_coeffs__')

        #lex = LineLex(open('deleteme.template', 'r'), '__standard_input_for_postprocess_coeffs_')
        lex.commenters = '#'            #(don't attempt to skip over comments)
        lex.line_extend_chars += '&'   #(because LAMMPS interprets '&' as '\')
        sys.stderr.write('word_terminators = "'+lex.word_terminators+'"\n')

        transition_rules_orig = []
        if_clause = []
        then_clause = []
        in_if_clause = True
        while lex:
            token = lex.get_token()
            if token == 'if':
                in_if_clause = True
                if ((len(if_clause)>0) and (len(then_clause)>0)):
                    transition_rules_orig.append([if_clause, then_clause])
                    if_clause = []
                    then_clause = []
                continue
            elif token == 'then':
                then_clause = []
                in_if_clause = False
                continue
            if in_if_clause:
                if_clause.append(token)
            else:
                then_clause.append(token)

        # now close the file (if we opened it)
        if not ('template' in args):
            templ_file.close()

        # Now split the if and then clauses into tokens separated by "and"
        if_conditions = []
        if_condition = []
        then_results = []
        then_result = []
        for rule in transition_rules_orig:
            if_clause = rule[0]
            then_clause = rule[1]
            #split the if_clause into lists of tokens separated by 'and':
            for token in if_clause:
                if ((token == 'and') and (len(if_condition) > 0)):
                    if_conditions.append(if_condition)
                else:
                    if_condition.append(token)
            if len(if_condition) > 0:
                if_conditions.append(if_condition)
            # Replace rule[0] with the if_conditions list
            rule[0] = if_conditions

            #split the then_clause into lists of tokens separated by 'and'
            for token in then_clause:
                if ((token == 'and') and (len(then_result) > 0)):
                    then_results.append(then_result)
                else:
                    then_result.append(token)
            if len(then_result) > 0:
                then_results.append(then_result)
            # Replace rule[1] with the then_results list
            rule[1] = then_results


        # Now loop through all of the transition rules.  For each rule,
        # figure out how many times the user specified an atom type, or
        # bonded type, or angle type or dihedral type or improper type
        # which must be satisfied in order to satisfy the if conditions.
        #
        # Then, for that rule, generate a separate transition rule for
        # every possible combination of atom types, bond types, angle types,
        # dihedral types, and improper types which satisfies the requirements
        # specified by the user after considering wildcards and regex.
        # In this way, a single rule specified by the user (with vague atom
        # type names or vague bonded type napes) might get translated
        # (expanded) into a large number of individual transition rules
        # for use with fix bond/react, where in each rule, each atom type,
        # bond type, angle type, etc... is specified explicitly.

        Natoms = 0 # we will store the number of atoms in the
                   # pre-reaction template here
        transition_rules = [] # we will store processed transition rules here
        atom_req = []  # atom type requirements
        bond_req = []  # bond type requirements
        angle_req = []  # angle type requirements
        dihedral_req = []  # dihedral type requirements
        improper_req = []  # improper type requirements
        for rule in transition_rules_orig:
            if_conditions = rule[0]
            for if_condition in if_conditions:
                tokens = if_condition
                assert(len(tokens) > 0)
                iatm = -1
                if tokens[0] == 'atom':
                    if not ((len(tokens) == 5) and
                            (tokens[1] == '[') and
                            (tokens[2].isnumeric() and
                             (int(tokens[2]) > 0)) and
                            (tokens[3] == ']') and
                            ((tokens[4].find('@') != -1) and
                             (tokens[4].find('atom:') != -1))):
                        raise InputError('Error in transitions file near:\n'+
                                         '      '+' '.join(tokens)+'\n')

                iatm = int(tokens[2])
                if iatm >= Natm:
                    atom_req += [] * (1 + iatm - Natm)
                    Natm = iatm + 1
                    assert(Natm == len(atom_req))

                typestr = tokens[4]  # a string identifying atom type(s)
                                     # eg: '@atom:/SPCE/O' or '@atom:C*'

                # If the second token is surrounded by '/' characters, interpret
                # it as a regular expression.
                type_is_re = HasRE(typestr)

                # If the second token contains wildcard characters, interpret
                # it as a wildcard (ie. glob) expression.
                type_is_wild = (HasWildcard(typestr) #does it contain *,?
                                and
                                (typestr[0] != '{')) #(ignore * or ? in {})

                if type_is_re:
                    regex_str = typestr[3:]
                    left_paren = text_after = ''
                    typepattern = re.compile(regex_str)
                else:
                    left_paren, typepattern, text_after = ExtractVarName(tokens[1])
                if (type_is_re or type_is_wild):
                    for atype in atom_types:
                        if MatchesPattern(atype, typepattern):
                            atom_req[iatm].append('@'+left_paren+atype+text_after)

        # ------------------ CONTINUEHERE --------------------
                    


    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()