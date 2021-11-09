#!/usr/bin/env python
import sys, math

g_filename    = __file__.split('/')[-1]
g_module_name  = g_filename
g_program_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str     = '2020-1-19'
g_version_str  = '0.3.0'



usage_example = g_program_name + """

merge_lines_periodic.py i1 i2 i3... \
          [-p natoms_per_monomer] \
          [-s nskip] [-d delim_atom] [-D delim_monomer] \
          < crds_input.raw  \
          > multiple_atom_crds_per_line.dat

Explanation: This script splits a text file into equally sized "blocks" (aka "monomers")
             and pastes the text text from different lines in each block into the 
             same line (with optional delimeters).
             The i1 i2 i3,... indices select the lines in each block (of atom
             coordinates in each monomer) that you want to merge together.
             Indexing begins at 0, not 1.  (The first line in a block has i=0)
         -Negative numbers correspond to atoms in the previous block(monomer).
         -Numbers larger than natoms_per_monomer lie in the next block(monomer).
             If any of these indices lie out of range, then the entire list 
             of lines in this block is ignored.
         -The -p argument indicates the number of lines in each block (aka "monomer")
             If the -p argument is skipped, then it is assumed to be infinity. (In other
             words, it is equal to the number of lines in the polymer conformation.)
         -The -s nskip argument allows you to skip over lines at the begining
             of the file.  (NOTE: Comments and lines beginning with comments 
             are ignored already, so don't include them in the nskip argument.)
         -The -d and -D delimeters allow you to change the string which 
             separates text belonging to different atoms(lines), and different 
             monomers (blocks). By default, they are " " and "\\n", respectively.
         -Blank lines (if present) in the input file are interpreted as delimeters 
             separating different "polymer conformations".  When encountered, each 
             "polymer conformation" is processed separately, with the output for 
             different polymer conformations delimted by blank lines.

"""



class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)


def EscCharStrToChar(s_in, escape='\\'):
    """ 
    EscCharStrToChar() replaces any escape sequences 
    in a string with their 1-character equivalents.

    """
    assert(len(escape) > 0)
    out_lstr = []
    escaped_state = False
    for c in s_in:
        if escaped_state:
            if (c == 'n'):
                out_lstr.append('\n')
            elif (c == 't'):
                out_lstr.append('\t')
            elif (c == 'r'):
                out_lstr.append('\r')
            elif (c == 'f'):
                out_lstr.append('\f')
            elif (c == '\''):
                out_lstr.append('\'')
            elif (c == '\"'):
                out_lstr.append('\"')
            elif c in escape:
                out_lstr.append(c)
            else:
                out_lstr.append(escape+c) # <- keep both characters
            escaped_state = False
        else:
            if c in escape:
                escaped_state = True
            else:
                out_lstr.append(c)

    return ''.join(out_lstr)



def SafelyEncodeString(in_str, 
                       quotes='\'\"',
                       delimiters=' \t\r\f\n', 
                       escape='\\', 
                       comment_char='#'):
    """
    SafelyEncodeString(in_str) scans through the input string (in_str),
    and returns a new string in which problematic characters 
    (like newlines, tabs, quotes, etc), are replaced by their two-character
    backslashed equivalents (like '\n', '\t', '\'', '\"', etc).  
    The escape character is the backslash by default, but it too can be
    overridden to create custom escape sequences 
    (but this does not effect the encoding for characters like '\n', '\t').

    """
    assert(len(escape) > 0)
    out_lstr = []
    use_outer_quotes = False
    for c in in_str:
        if (c == '\n'):
            c = '\\n'
        elif (c == '\t'):
            c = '\\t'
        elif (c == '\r'):
            c = '\\r'
        elif (c == '\f'):
            c = '\\f'
        elif c in quotes:
            c = escape[0]+c
        elif c in escape:
            c = c+c
        elif c in delimiters:
            use_outer_quotes = True
        # hmm... that's all that comes to mind.  Did I leave anything out?
        out_lstr.append(c)

    if use_outer_quotes:
        out_lstr = ['\"'] + out_lstr + ['\"']

    return ''.join(out_lstr)



def ProcessSnapshot(lines, 
                    out_file, 
                    offsets, 
                    period, 
                    nskip, 
                    delimeter_atom,
                    delimeter_monomer):

    offsets_min = min(offsets)
    offsets_max = max(offsets)
    if period == 0:
        num_monomers = 1
    else:
        num_monomers = math.floor((len(lines)-nskip)/period)
    for I in range(0, num_monomers):

        # If any of the entries will be missing, then ignore the whole list
        # of atoms (lines) for this monomer (block).
        if (I*period + offsets_min < nskip):
            continue
        if (I*period + offsets_max >= len(lines)):
            continue

        for J in range(0, len(offsets)):
            j = offsets[J]
            i = (I*period + nskip) + j
            if (nskip <= i) and (i < len(lines)):
                out_file.write(lines[i])
                if J+1 < len(offsets):
                    out_file.write(delimeter_atom)
                else:
                    out_file.write(delimeter_monomer)


g_period = 0
g_nskip = 0
g_delimeter_atom = ' '
g_delimeter_monomer = '\n'
g_delimeter_snapshot = '\n'
g_offsets = []


#######  Main Code Below: #######

sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
if sys.version < '3':
    sys.stderr.write(' (python version < 3)\n')
else:
    sys.stderr.write('\n')

try:

    argv = [arg for arg in sys.argv]

    # Loop over the remaining arguments not processed yet.
    # These arguments are specific to the lttree.py program
    # and are not understood by ttree.py:
    i = 1
    while i < len(argv):
        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
        if ((argv[i].lower() == '-?') or
            (argv[i].lower() == '--?') or
            (argv[i].lower() == '-help') or
            (argv[i].lower() == '-help')):
            if i+1 >= len(argv):
                sys.stdout.write("\n Usage:\n\n"+usage_example+'\n')
                sys.exit(0)

        elif argv[i].lower() == '-p':
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a number.\n')
            g_period = int(argv[i+1])
            sys.stderr.write('   period = '+str(g_period)+'\n')
            del(argv[i:i+2])

        elif argv[i].lower() == '-s':
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a number.\n')
            g_nskip = float(argv[i+1])
            sys.stderr.write('   skip first '+str(g_nskip)+' non-comment lines\n')
            del(argv[i:i+2])

        elif argv[i].lower() == '-d':
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a string.\n')
            g_delimeter_atom = EscCharStrToChar(argv[i+1])
            sys.stderr.write('   delimeter_atom = \"'+SafelyEncodeString(g_delimeter_atom)+'\"\n')
            del(argv[i:i+2])

        elif argv[i].lower() == '-D':
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by string.\n')
            g_delimeter_atom = EscCharStrToChar(argv[i+1])
            sys.stderr.write('   delimeter_monomer = \"'+SafelyEncodeString(g_delimeter_monomer)+'\"\n')
            del(argv[i:i+2])

        elif argv[i][0] == '-':
            # Note: It could be a negative integer, so check for 
            # that before printing an error message
            if not argv[i][1:].isdigit():
                raise InputError('Error('+g_program_name+'):\n'
                                 'Unrecogized command line argument \"'+argv[i]+'\"\n')
            i += 1

        else:
            i += 1

    if len(argv) == 1:
        raise InputError("Error: Expected a list of integers.\n\n"+
                         "Usage: \n\n"+
                         "       "+usage_example+"\n")

    g_offsets = [int(arg) for arg in argv[1:]]





    # --- Now (finally) read the lines in the standard input ----
    n_snapshots = 0
    lines = []
    in_file = sys.stdin
    for line_orig in in_file:

        ic = line_orig.find('#')
        if ic != -1:
            line = line_orig[:ic]
        else:
            line = line_orig.rstrip('\n')    

        # Blank lines in a trajectory file usually signal the end of the
        # coordinates for that snapshot in the trajectory, and the beginning
        # of the next snapshot.

        if len(line_orig.strip()) == 0:
            if n_snapshots > 0:
                sys.stdout.write(g_delimeter_snapshot)
            if len(lines) > 0:
                ProcessSnapshot(lines, 
                                sys.stdout,
                                g_offsets, 
                                g_period,
                                g_nskip, 
                                g_delimeter_atom,
                                g_delimeter_monomer)
                n_snapshots += 1
            # Clear the lines buffer to begin reading the new snapshot
            del lines[:]
        else:
            if len(line.strip()) > 0:
                lines.append(line)

    if len(lines) > 0:
        if n_snapshots > 0:
            sys.stdout.write(g_delimeter_snapshot)
        # After reading all of the lines in the file, deal with any lines 
        # left over since reading the last frame
        ProcessSnapshot(lines, 
                        sys.stdout,
                        g_offsets, 
                        g_period,
                        g_nskip, 
                        g_delimeter_atom,
                        g_delimeter_monomer)







except (ValueError, InputError) as err:
    sys.stderr.write('\n\n'+str(err)+'\n')
    sys.exit(-1)



