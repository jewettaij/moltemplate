# -*- coding: iso-8859-1 -*-
### -*- coding: utf-8 -*-

# Authors: Eric S. Raymond, 21 Dec 1998
#          Andrew Jewett (jewett.aij at g mail)
# LICENSE: The PSF license:
# https://docs.python.org/3/license.html
# The PSF license is compatible with the GPL license.  It is not a copyleft
# license.  It is apparently similar to the BSD and MIT licenses.
#
# Contributions:
# Module and documentation by Eric S. Raymond, 21 Dec 1998
# Input stacking and error message cleanup added by ESR, March 2000
# push_source() and pop_source() made explicit by ESR, January 2001.
# Posix compliance, split(), string arguments, and
# iterator interface by Gustavo Niemeyer, April 2003.
# Unicode support hack ("wordterminators") and numerous other hideous
# ttree-specific hacks added by Andrew Jewett September 2011.


"""A lexical analyzer class for simple shell-like syntaxes.
   This version has been modified slightly to work better with unicode.
   It was forked from the version of shlex that ships with python 3.2.2.
   A few minor features and functions have been added.  -Andrew Jewett 2011 """


import os.path
import sys
from collections import deque
import re
import fnmatch
import string
#import gc


try:
    from cStringIO import StringIO
except ImportError:
    try:
        from StringIO import StringIO
    except ImportError:
        from io import StringIO

__all__ = ["TtreeShlex",
           "split",
           "LineLex",
           "SplitQuotedString",
           "ExtractVarName",
           "GetVarName",
           "EscCharStrToChar",
           "SafelyEncodeString",
           "RemoveOuterQuotes",
           "MaxLenStr",
           "VarNameToRegex",
           "HasRE",
           "HasWildcard",
           "MatchesPattern",
           "InputError",
           "ErrorLeader",
           "SrcLoc",
           "OSrcLoc",
           "TextBlock",
           "VarRef",
           "VarNPtr",
           "VarBinding",
           "SplitTemplate",
           "SplitTemplateMulti",
           "TableFromTemplate",
           "ExtractCatName",
           #"_TableFromTemplate",
           #"_DeleteLineFromTemplate",
           "DeleteLinesWithBadVars",
           "TemplateLexer"]


class TtreeShlex(object):
    """ A lexical analyzer class for simple shell-like syntaxes.
    TtreeShlex is a backwards-compatible version of python's standard shlex
    module. It has the additional member: "self.wordterminators", which
    overrides the "self.wordchars" member.  This enables better handling of
    unicode characters by allowing a much larger variety of characters to
    appear in words or tokens parsed by TtreeShlex.

    """

    def __init__(self,
                 instream=None,
                 infile=None,
                 posix=False):
        if isinstance(instream, str):
            instream = StringIO(instream)
        if instream is not None:
            self.instream = instream
            self.infile = infile
        else:
            self.instream = sys.stdin
            self.infile = None
        self.posix = posix
        if posix:
            self.eof = None
        else:
            self.eof = ''
        self.commenters = '#'
        self.wordchars = ('abcdfeghijklmnopqrstuvwxyz'
                          'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')
        #if self.posix:
        #    self.wordchars += ('��������������������������������'
        #                       '������������������������������')

        if self.posix:
            self.wordchars += ('��������������������������������'
                               '������������������������������')

        self.wordterminators = set([])
        self.prev_space_terminator = ''
        self.whitespace = ' \t\r\f\n'
        self.whitespace_split = False
        self.quotes = '\'"'
        self.escape = '\\'
        self.escapedquotes = '"'
        self.operators = '='  #binary numeric operators like +-*/ might be added
        self.state = ' '
        self.pushback = deque()
        self.lineno = 1
        self.debug = 0
        self.token = ''
        self.filestack = deque()
        # self.source_triggers
        # are tokens which allow the seamless insertion of other
        # files into the file being read.
        self.source_triggers = set(['source'])
        self.source_triggers_x = set([])
        # self.source_triggers_x is a subset of self.source_triggers.
        # In this case file inclusion is exclusive.
        # In other words, the file is only included
        # if it has not been included already.  It does this
        # by checking if one of these tokens has been encountered.
        self.source_files_restricted = set([])
        self.include_path = []
        if 'TTREE_PATH' in os.environ:
            include_path_list = os.environ['TTREE_PATH'].split(':')
            self.include_path += [d for d in include_path_list if len(d) > 0]
        if self.debug:
            sys.stderr.write('TtreeShlex: reading from %s, line %d'
                  % (self.instream, self.lineno))
        self.end_encountered = False

    @staticmethod
    def _belongs_to(char, include_chars, exclude_chars):
        if ((not exclude_chars) or (len(exclude_chars)==0)):
            return char in include_chars
        else:
            return char not in exclude_chars

    def push_raw_text(self, text):
        """Push a block of text onto the stack popped by the ReadLine() method.
        (If multiple lines are present in the text, (which is determined by
        self.line_terminators), then the text is split into multiple lines
        and each one of them is pushed onto this stack individually.
        The "self.lineno" counter is also adjusted, depending on the number
        of newline characters in "line".
            Do not strip off the newline, or other line terminators
            at the end of the text block before using push_raw_text()!

        """
        if self.debug >= 1:
            sys.stderr.write("TtreeShlex: pushing token " + repr(text))
        for c in reversed(text):
            self.pushback.appendleft(c)
            if c == '\n':
                self.lineno -= 1
        if len(text) > 0:
            self.end_encountered = False

    def push_token(self, text):
        "Push a token onto the stack popped by the get_token method"
        self.push_raw_text(text + self.prev_space_terminator)

    def push_source(self, newstream, newfile=None):
        "Push an input source onto the lexer's input source stack."
        if isinstance(newstream, str):
            newstream = StringIO(newstream)
        self.filestack.appendleft((self.infile, self.instream, self.lineno))
        self.infile = newfile
        self.instream = newstream
        self.lineno = 1
        if self.debug:
            if newfile is not None:
                sys.stderr.write('TtreeShlex: pushing to file %s' % (self.infile,))
            else:
                sys.stderr.write('TtreeShlex: pushing to stream %s' % (self.instream,))

    def pop_source(self):
        "Pop the input source stack."
        self.instream.close()
        (self.infile, self.instream, self.lineno) = self.filestack.popleft()
        if self.debug:
            sys.stderr.write('TtreeShlex: popping to %s, line %d'
                  % (self.instream, self.lineno))
        self.state = ' '

    def get_token(self):
        "Get a token from the input stream (or from stack if it's nonempty)"
        #### #CHANGING: self.pushback is now a stack of characters, not tokens
        #### if self.pushback:
        ####    tok = self.pushback.popleft()
        ####    if self.debug >= 1:
        ####        sys.stderr.write("TtreeShlex: popping token " + repr(tok))
        ####    return tok
        #### No pushback.  Get a token.
        raw = self.read_token()
        # Handle inclusions
        if self.source_triggers is not None:
            while raw in self.source_triggers:
                fname = self.read_token()
                spec = self.sourcehook(fname)
                if spec:
                    (newfile, newstream) = spec
                    if ((raw not in self.source_triggers_x) or
                            (newfile not in self.source_files_restricted)):
                        self.push_source(newstream, newfile)
                        if raw in self.source_triggers_x:
                            self.source_files_restricted.add(newfile)
                    else:
                        if self.debug >= 1:
                            sys.stderr.write(
                                '\ndebug warning: duplicate attempt to import file:\n               \"' + newfile + '\"\n')
                raw = self.get_token()

        # Maybe we got EOF instead?
        while raw == self.eof:
            if not self.filestack:
                return self.eof
            else:
                self.pop_source()
                raw = self.get_token()
        # Neither inclusion nor EOF
        if self.debug >= 1:
            if raw != self.eof:
                sys.stderr.write("TtreeShlex: token=" + repr(raw))
            else:
                sys.stderr.write("TtreeShlex: token=EOF")

        if raw == self.eof:
            self.end_encountered = True

        return raw

    def read_char(self):
        if self.pushback:
            nextchar = self.pushback.popleft()
            assert((type(nextchar) is str) and (len(nextchar)==1))
        else:
            nextchar = self.instream.read(1)
        return nextchar

    def read_token(self):
        self.prev_space_terminator = ''
        quoted = False
        escapedstate = ' '
        while True:
            #### self.pushback is now a stack of characters, not tokens
            nextchar = self.read_char()
            if nextchar == '\n':
                self.lineno = self.lineno + 1
            if self.debug >= 3:
                sys.stderr.write("TtreeShlex: in state", repr(self.state),
                      "I see character:", repr(nextchar))
            if self.state is None:
                self.token = ''        # past end of file
                break
            elif self.state == ' ':
                if not nextchar:
                    self.state = None  # end of file
                    break
                elif nextchar in self.whitespace:
                    if self.debug >= 2:
                        sys.stderr.write("TtreeShlex: I see whitespace in whitespace state")
                    if self.token or (self.posix and quoted):
                        # Keep track of which whitespace
                        # character terminated the token.
                        self.prev_space_terminator = nextchar
                        break   # emit current token
                    else:
                        continue
                elif nextchar in self.commenters:
                    self.instream.readline()
                    self.lineno = self.lineno + 1
                elif self.posix and nextchar in self.escape:
                    escapedstate = 'a'
                    self.state = nextchar
                elif TtreeShlex._belongs_to(nextchar,
                                            self.wordchars,
                                            self.wordterminators):
                    self.token = nextchar
                    self.state = 'a'
                elif nextchar in self.quotes:
                    if not self.posix:
                        self.token = nextchar
                    self.state = nextchar
                elif self.whitespace_split:
                    self.token = nextchar
                    self.state = 'a'
                else:
                    self.token = nextchar
                    if self.token or (self.posix and quoted):
                        break   # emit current token
                    else:
                        continue
            elif self.state in self.quotes:
                quoted = True
                if not nextchar:      # end of file
                    if self.debug >= 2:
                        sys.stderr.write("TtreeShlex: I see EOF in quotes state")
                    # XXX what error should be raised here?
                    raise ValueError("Error at or before " + self.error_leader() + "\n"
                                     "      No closing quotation.")
                if nextchar == self.state:
                    if not self.posix:
                        self.token = self.token + nextchar
                        self.state = ' '
                        break
                    else:
                        self.state = 'a'
                elif self.posix and nextchar in self.escape and \
                        self.state in self.escapedquotes:
                    escapedstate = self.state
                    self.state = nextchar
                else:
                    self.token = self.token + nextchar
            elif self.state in self.escape:
                if not nextchar:      # end of file
                    if self.debug >= 2:
                        sys.stderr.write("TtreeShlex: I see EOF in escape state")
                    # What error should be raised here?
                    raise InputError('File terminated immediately following an escape character.')
                # In posix shells, only the quote itself or the escape
                # character may be escaped within quotes.
                if escapedstate in self.quotes and \
                   nextchar != self.state and nextchar != escapedstate:
                    self.token = self.token + self.state
                self.token = self.token + nextchar
                self.state = escapedstate
            elif self.state == 'a':
                if not nextchar:
                    self.state = None   # end of file
                    break
                elif nextchar in self.whitespace:
                    if self.debug >= 2:
                        sys.stderr.write("TtreeShlex: I see whitespace in word state")
                    self.state = ' '
                    if self.token or (self.posix and quoted):
                        # Keep track of which whitespace
                        # character terminated the token.
                        self.prev_space_terminator = nextchar
                        break   # emit current token
                    else:
                        continue
                elif nextchar in self.commenters:
                    comment_contents = self.instream.readline()
                    self.lineno = self.lineno + 1
                    if self.posix:
                        self.state = ' '
                        if self.token or (self.posix and quoted):
                            # Keep track of which character(s) terminated
                            # the token (including whitespace and comments).
                            self.prev_space_terminator = nextchar + comment_contents
                            break   # emit current token
                        else:
                            continue
                elif self.posix and nextchar in self.quotes:
                    self.state = nextchar
                elif self.posix and nextchar in self.escape:
                    escapedstate = 'a'
                    self.state = nextchar
                elif (TtreeShlex._belongs_to(nextchar,
                                             self.wordchars,
                                             self.wordterminators)
                      or (nextchar in self.quotes)
                      or (self.whitespace_split)):
                    self.token = self.token + nextchar
                else:
                    self.pushback.appendleft(nextchar)
                    if self.debug >= 2:
                        sys.stderr.write("TtreeShlex: I see punctuation in word state")
                    self.state = ' '
                    if self.token:
                        break   # emit current token
                    else:
                        continue
        result = self.token
        self.token = ''
        if self.posix and not quoted and result == '':
            result = None
        if self.debug > 1:
            if result:
                sys.stderr.write("TtreeShlex: raw token=" + repr(result))
            else:
                sys.stderr.write("TtreeShlex: raw token=EOF")
        return result

    def sourcehook(self, newfile):
        "Hook called on a filename to be sourced."
        newfile = RemoveOuterQuotes(newfile)
        # This implements cpp-like semantics for relative-path inclusion.
        newfile_full = newfile
        if isinstance(self.infile, str) and not os.path.isabs(newfile):
            newfile_full = os.path.join(os.path.dirname(self.infile), newfile)
        try:
            f = open(newfile_full, "r")
        except IOError:
            # If not found,
            err = True
            # ...then check to see if the file is in one of the
            # directories in the self.include_path list.
            for d in self.include_path:
                newfile_full = os.path.join(d, newfile)
                try:
                    f = open(newfile_full, "r")
                    err = False
                    break
                except IOError:
                    err = True
            if err:
                raise InputError('Error at ' + self.error_leader() + '\n'
                                 '       unable to open file \"' + newfile + '\"\n'
                                 '       for reading.\n')
        return (newfile, f)

    def error_leader(self, infile=None, lineno=None):
        "Emit a C-compiler-like, Emacs-friendly error-message leader."
        if infile is None:
            infile = self.infile
        if lineno is None:
            lineno = self.lineno
        return "\"%s\", line %d: " % (infile, lineno)

    def __iter__(self):
        return self

    def __next__(self):
        token = self.get_token()
        if token == self.eof:
            raise StopIteration
        return token

    def __bool__(self):
        return not self.end_encountered

    # For compatibility with python 2.x, I must also define:
    def __nonzero__(self):
        return self.__bool__()


# The split() function was originally from shlex
# It is included for backwards compatibility.
def split(s, comments=False, posix=True):
    lex = TtreeShlex(s, posix=posix)
    lex.whitespace_split = True
    if not comments:
        lex.commenters = ''
    return list(lex)


##################### NEW ADDITIONS (may be removed later) #################

#"""
#  -- linelex.py --
# linelex.py defines the LineLex class, which inherits from, and further
# augments the capabilities of TtreeShlex by making it easier to parse
# individual lines one at a time.  (The original shlex's "source" inclusion
# ability still works when reading entire lines, and lines are still counted.)
#
#"""

#import sys


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


def ErrorLeader(infile, lineno):
    return '\"' + infile + '\", line ' + str(lineno)


class SrcLoc(object):
    """ SrcLoc is essentially nothing more than a 2-tuple containing the name
    of a file (str) and a particular line number inside that file (an integer).

    """
    __slots__ = ["infile", "lineno"]

    def __init__(self, infile='', lineno=-1):
        self.infile = infile
        self.lineno = lineno


def SplitQuotedString(string,
                      quotes='\'\"',
                      delimiters=' \t\r\f\n',
                      escape='\\',
                      comment_char='#',
                      endquote=None):
    tokens = []
    token = ''
    reading_token = True
    escaped_state = False
    quote_state = None
    for c in string:

        if (c in comment_char) and (not escaped_state) and (quote_state == None):
            tokens.append(token)
            return tokens

        elif (c in delimiters) and (not escaped_state) and (quote_state == None):
            if reading_token:
                tokens.append(token)
                token = ''
                reading_token = False

        elif c in escape:
            if escaped_state:
                token += c
                reading_token = True
                escaped_state = False
            else:
                escaped_state = True
                # and leave c (the '\' character) out of token
        elif (c == quote_state) and (not escaped_state) and (quote_state != None):
            quote_state = None
            if include_endquote:
                token += c
        elif (c in quotes) and (not escaped_state):
            if quote_state == None:
                if endquote != None:
                    quote_state = endquote
                else:
                    quote_state = c
                # Now deal with strings like
                #    a "b" "c d" efg"h i j"
                # Assuming quotes='"', then we want this to be split into:
                #    ['a', 'b', 'c d', 'efg"h i j"']
                # ...in other words, include the end quote if the token did
                #    not begin with a quote
                include_endquote = False
                if token != '':
                    # if this is not the first character in the token
                    include_endquote = True
            token += c
            reading_token = True
        else:
            if (c == 'n') and (escaped_state == True):
                c = '\n'
            elif (c == 't') and (escaped_state == True):
                c = '\t'
            elif (c == 'r') and (escaped_state == True):
                c = '\r'
            elif (c == 'f') and (escaped_state == True):
                c = '\f'
            token += c
            reading_token = True
            escaped_state = False

    # Remove any empty strings from the front or back of the list,
    # just in case SplitQuotedString() fails to remove them.
    # (Possible bug in SplitQuotedString(), but too lazy to investigate.)
    if (len(tokens) > 0) and (tokens[0] == ''):
        del tokens[0]
    if (len(tokens) > 0) and (tokens[-1] == ''):
        del tokens[-1]

    if (len(string) > 0) and (token != ''):
        tokens.append(token)
    return tokens




def GetVarName(lex):
    """ Read a string like 'atom:A  '  or  '{/atom:A B/C/../D }ABC '
        and return ('','atom:A','  ')  or  ('{','/atom:A B/C/../D ','}ABC')
        These are 3-tuples containing the portion of the text containing 
        only the variable's name (assumed to be within the text),
        ...in addition to the text on either side of the variable name.
    """
    escape = '\''
    lparen = '{'
    rparen = '}'
    if hasattr(lex, 'escape'):
        escape = lex.escape
    if hasattr(lex, 'var_open_paren'):
        lparen = lex.var_open_paren
    if hasattr(lex, 'var_close_paren'):
        rparen = lex.var_close_paren

    nextchar = lex.read_char()
    # Skip past the left-hand side paren '{'
    paren_depth = 0
    escaped = False
    if nextchar == lparen:
        paren_depth = 1
    elif nextchar in lex.escape:
        escaped = True
    elif (hasattr(lex, 'wordterminators') and
          (nextchar in lex.wordterminators)):
        lex.push_raw_text(nextchar)
        return ''
    else:
        lex.push_raw_text(nextchar)
    # Now read the variable name:
    var_name_l = []
    while lex:
        nextchar=lex.read_char()
        if nextchar == '':
            break
        elif nextchar == '\n':
            lex.lineno += 1
            if paren_depth > 0:
                var_name_l.append(nextchar)
            else:
                lex.push_raw_text(nextchar)
                break
        elif escaped:
            var_name_l.append(nextchar)
            escaped = False
        elif nextchar in lex.escape:
            escaped = True
        elif nextchar == lparen:
            paren_depth += 1
            if (hasattr(lex, 'wordterminators') and
                (nextchar in lex.wordterminators)):
                lex.push_raw_text(nextchar)
                break
            else:
                var_name_l.append(nextchar)
        elif nextchar == rparen:
            paren_depth -= 1
            if paren_depth == 0:
                break
            elif (hasattr(lex, 'wordterminators') and
                  (nextchar in lex.wordterminators)):
                lex.push_raw_text(nextchar)
                break
            else:
                var_name_l.append(nextchar)
        elif paren_depth > 0:
            var_name_l.append(nextchar)
            escaped = False
        elif nextchar in lex.whitespace:
            lex.push_raw_text(nextchar)
            break
        elif (hasattr(lex, 'wordterminators') and
              (nextchar in lex.wordterminators) and
              (paren_depth == 0)):
            lex.push_raw_text(nextchar)
            break
        elif nextchar in lex.commenters:
            lex.instream.readline()
            lex.lineno += 1
            break
        else:
            var_name_l.append(nextchar)
            escaped = False
    var_name = ''.join(var_name_l)
    return var_name



def ExtractVarName(text,
                   commenters = '#',
                   whitespace = ' \t\r\f\n'):
    """ Read a string like 'atom:A  '  or  '{/atom:A B/C/../D }ABC '
        and return ('','atom:A','  ')  or  ('{','/atom:A B/C/../D ','}ABC')
        These are 3-tuples containing the portion of the text containing 
        only the variable's name (assumed to be within the text),
        ...in addition to the text on either side of the variable name.
    """
    ibegin = 0
    left_paren = ''
    if text[0] == '{':
        ibegin = 1
        left_paren = text[0] #(GetVarName() strips the leading '{' character)
    # The best way to insure consistency with other code is to use
    # lex.GetVarName() to figure out where the variable name ends.
    lex = TtreeShlex(StringIO(text))
    var_name = GetVarName(lex)
    # Any text following the end of the variable name should be returned as well
    text_after_list = []
    if left_paren:
        text_after_list.append('}') #(GetVarName() strips the trailing '}' char)
    while lex:
        c = lex.read_char()
        if c == '':
            break
        text_after_list.append(c)
    text_after = ''.join(text_after_list)
    return (left_paren, var_name, text_after)


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
                out_lstr.append(escape + c)  # <- keep both characters
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
    and returns a new string in which probletic characters
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
            c = escape[0] + c
        elif c in escape:
            c = c + c
        elif c in delimiters:
            use_outer_quotes = True
        # hmm... that's all that comes to mind.  Did I leave anything out?
        out_lstr.append(c)

    if use_outer_quotes:
        out_lstr = ['\"'] + out_lstr + ['\"']

    return ''.join(out_lstr)


def RemoveOuterQuotes(text, quotes='\"\''):
    if ((len(text) >= 2) and (text[0] in quotes) and (text[-1] == text[0])):
        return text[1:-1]
    else:
        return text


def MaxLenStr(s1, s2):
    if len(s2) > len(s1):
        return s2
    else:
        return s1


def VarNameToRegex(s):
    """
    Returns the portion of a TTREE-style variable name (eg "@atom:re.C[1-5]")
    that corresponds to a regular expression (eg "C[1-5]").  A variable name
    is assumed to encode a regular expression if it begins with "re.", OR if
    the a ':' character is followed by "re.".
    If so, the text in s (excluding "re.") is assumed to be a regular expresion
    and is returned to the caller.
    If not, the empty string ('') is returned.
    If the first or second character is a '{', and if the final character
    is '}', they will be deleted.  Consequently:
      VarNameToRegex('@atom:C') returns ''
      VarNameToRegex('@atom:re.C[1-5]') returns '@atom:C[1-5]'
      VarNameToRegex('@{/atom:re.C[1-5]}') returns '@/atom:C[1-5]'
      VarNameToRegex('@bond:AB') returns ''
      VarNameToRegex('@bond:re.A*B') returns '@bond:a*b'
      VarNameToRegex('bond:re.A*B') returns 'bond:a*b'
      VarNameToRegex('{bond:re.A*B}') returns 'bond:a*b'
      VarNameToRegex('@{bond:re.A*B}') returns '@bond:a*b'
    """
    # First, deal with parenthesis {}
    iparen_L = s.find('{')
    iparen_R = s.rfind('}')
    if (((iparen_L == 0) or (iparen_L == 1)) and (iparen_R == len(s)-1)):
        optional_char = ''
        if iparen_L == 1:
            optional_char = s[0]
        s = optional_char + s[iparen_L+1:iparen_R]
    # Now check to see if the remaining string contains 're.' or ':re.'
    icolon = s.find(':')
    # If 're.' is not found immediately after the first ':' character
    # or following a '/' character
    # (or if it is not found at the beginning when no ':' is present)
    # then there is no regular expression.  In that case, return ''
    ire = s.find('re.')
    if ((ire == -1) or
        (not ((ire > 0) and ((s[ire-1] == ':') or (s[ire-1] == '/'))))):
        return ''
    return s[0:ire] + s[ire+3:]
    

def HasRE(pat):
    """
    Returns true if a string (pat) begins with 're.'
    """
    return len(VarNameToRegex(pat)) > 0


def HasWildcard(pat):
    """
    Returns true if a string (pat) contains a '*' or '?' character.

    """
    return (pat.find('*') != -1) or (pat.find('?') != -1)


# def HasWildcard(pat):
#    """
#    Returns true if a string (pat) contains a non-backslash-protected
#    * or ? character.
#
#    """
#    N=len(pat)
#    i=0
#    while i < N:
#        i = pat.find('*', i, N)
#        if i == -1:
#            break
#        elif (i==0) or (pat[i-1] != '\\'):
#            return True
#        i += 1
#    i=0
#    while i < N:
#        i = pat.find('?', i, N)
#        if i == -1:
#            break
#        elif (i==0) or (pat[i-1] != '\\'):
#            return True
#        i += 1
#    return False


def MatchesPattern(s, pattern):
    if type(pattern) is str:
        # old code:
        # if ((len(s) > 1) and (s[0] == '/') and (s[-1] == '/'):
        #    re_string = p[1:-1]  # strip off the slashes '/' and '/'
        #    if not re.search(re_string, s):
        #        return False
        # new code:
        #    uses precompiled regular expressions (See "pattern.search" below)
        if HasWildcard(pattern):
            if not fnmatch.fnmatchcase(s, pattern):
                return False
        elif s != pattern:
            return False
    else:
        #assert(type(p) is _sre.SRE_Match)
        # I assume pattern = re.compile(some_reg_expr)
        if not pattern.search(s):
            return False
    return True


def MatchesAll(multi_string, pattern):
    assert(len(multi_string) == len(pattern))
    for i in range(0, len(pattern)):
        if not MatchesPattern(multi_string[i], pattern[i]):
            return False
    return True


class LineLex(TtreeShlex):
    """ This class extends the TtreeShlex module (a slightly modified
    version of the python 3.2.2 version of shlex).  LineLex has the
    ability to read one line at a time (in addition to one token at a time).
    (Many files and scripts must be parsed one line at a time instead of one
     token at a time.  In these cases, the whitespace position also matters.)

    Arguably, this class might not be necessary.
    I could get rid of this class completely.  That would be nice.  To do that
    we would need to augment and generalize shlex's get_token() member function
    to make it read lines, not just tokens.  Of course, you can always
    change the wordchars (or wordterminators).  Even so, there are two other
    difficulties using the current version of shlex.get_token() to read lines:
    1) File inclusion happen whenever the beginning of a line/token matches one
       of the "source_triggers" (not the whole line as required by get_token()).
    2) Lines ending in a special character (by default the backslash character)
       continue on to the next line.
    This code seems to work on our test files, but I'm sure there are bugs.
    Andrew 2012-3-25

    """

    def __init__(self,
                 instream=None,
                 infile=None,
                 posix=False):
        TtreeShlex.__init__(self, instream, infile, posix)
        self.line_terminators = '\n'
        self.line_extend_chars = '\\'
        self.skip_comments_during_readline = True

    def _StripComments(self, line):
        if self.skip_comments_during_readline:
            for i in range(0, len(line)):
                if ((line[i] in self.commenters) and
                        ((i == 0) or (line[i - 1] not in self.escape))):
                    return line[:i]
        return line

    def _ReadLine(self,
                  recur_level=0):
        """
        This function retrieves a block of text, halting at a
        terminal character.  Escape sequences are respected.
        The self.lineno (newline counter) is also maintained.

        The main difference between Readline and get_token()
        is the way they handle the "self.source_triggers" member.
        Both Readline() and get_token() insert text from other files when they
        encounter a string in "self.source_triggers" in the text they read.
        However ReadLine() ONLY inserts text from other files if the token which
        matches with self.source_triggers appears at the beginning of the line.
        get_token() inserts text only if lex.source matches the entire token.

        comment-to-self:
         At some point, once I'm sure this code is working, I should replace
         shlex.get_token() with the code from ReadLine() which is more general.
         It would be nice to get rid of "class LineLex" entirely.  ReadLine()
         is the only new feature that LineLex which was lacking in shlex.

         To do this I would need to add a couple optional arguments to
         "get_token()", allowing it to mimic ReadLine(), such as:
           "override_wordterms" argument (which we can pass a '\n'), and
           "token_extender" argument (like '\' for extending lines)

        """
        first_token = ''
        line = ''
        escaped_state = False
        found_space = False
        while True:
            nextchar = self.read_char()
            # sys.stderr.write('nextchar=\"'+nextchar+'\"\n')
            while nextchar == '':
                if not self.filestack:
                    return self._StripComments(line), '', first_token, found_space
                else:
                    self.pop_source()
                    nextchar = self.read_char()
            if nextchar == '\n':
                self.lineno += 1

            if escaped_state:
                escaped_state = False
            else:
                if nextchar in self.escape:
                    line += nextchar
                    escaped_state = True
                else:
                    escaped_state = False

            if not escaped_state:
                if (nextchar in self.whitespace):
                    found_space = True
                    while first_token in self.source_triggers:
                        fname = RemoveOuterQuotes(self.get_token())
                        if (fname == '') or (fname in self.source_triggers):
                            raise InputError('Error: near ' + self.error_leader() + '\n'
                                             '       Nonsensical file inclusion request.\n')
                        if self.debug >= 0:
                            sys.stderr.write(('  ' * recur_level) +
                                             'reading file \"' + fname + '\"\n')
                        spec = self.sourcehook(fname)
                        if spec:
                            (fname, subfile) = spec
                            if ((first_token not in self.source_triggers_x) or
                                (fname not in self.source_files_restricted)):
                                self.push_source(subfile, fname)
                            if first_token in self.source_triggers_x:
                                self.source_files_restricted.add(fname)
                            else:
                                if self.debug >= 0:
                                    sys.stderr.write('\nWarning at ' + self.error_leader() + ':\n'
                                                     '          duplicate attempt to import file:\n         \"' + fname + '\"\n')

                        line, nextchar, first_token, found_space = \
                            self._ReadLine(recur_level + 1)

                if nextchar in self.line_terminators:
                    line_nrw = line.rstrip(self.whitespace)
                    # sys.stderr.write('line_nrw=\"'+line_nrw+'\"\n')
                    if ((len(line_nrw) > 0) and
                            (line_nrw[-1] in self.line_extend_chars) and
                            ((len(line_nrw) < 2) or (line_nrw[-2] not in self.escape))):
                        # delete the line_extend character
                        line = line_nrw[:-1]
                        # from the end of that line and keep reading...
                    else:
                        return self._StripComments(line), nextchar, first_token, found_space
                else:
                    line += nextchar
                    if not found_space:
                        first_token += nextchar

    def ReadLine(self, recur_level=0):
        line, nextchar, first_token, found_space = \
            self._ReadLine(recur_level)
        if nextchar == self.eof:
            self.end_encountered = True
        return line + nextchar

    @staticmethod
    def TextBlock2Lines(text, delimiters, keep_delim=True):
        """ This splits a string into a list of sub-strings split by delimiter
        characters.  This function is different from the standard str.split()
        function: The string is split at every character which belongs to the
        "delimiters" argument (which can be a string or some other container).
        This character is included at the end of every substring.  Example:
        TextBlock2Lines('\nabc\nde^fg\nhi j\n', '^\n')
        returns:
        ['\n', 'abc\n', 'de^', 'fg\n', 'hi j\n']

        """
        ls = []
        i = 0
        i_prev = 0
        while i < len(text):
            if text[i] in delimiters:
                if keep_delim:
                    ls.append(text[i_prev:i + 1])
                else:
                    ls.append(text[i_prev:i])
                i_prev = i + 1
            i += 1
        if (i_prev < len(text)):
            ls.append(text[i_prev:i + 1])
        return ls

    def __iter__(self):
        return self

    def __next__(self):
        line = self.ReadLine()
        if line == self.eof:
            raise StopIteration
        return line


class OSrcLoc(object):
    """ OSrcLoc is barely more than a 2-tuple containing the name of a file
        (a string) and a particular line number inside that file (an integer).
        These objects are passed around and stored in the nodes of
        every tree, so that if a syntax error or broken link in that node
        is discovered, an error message can be provided to the user.

    """

    __slots__ = ["infile", "lineno", "order"]
    count = 0

    def __init__(self, infile='', lineno=-1):
        self.infile = infile
        self.lineno = lineno
        OSrcLoc.count += 1
        self.order = OSrcLoc.count  # keep track of how many times it was called

    def __lt__(self, x):
        return self.order < x.order

    # def __repr__(self):
    #    return repr((self.infile, self.lineno, self.order))


class TextBlock(object):
    """TextBlock is just a 3-tuple consisting of a string, and an OSrcLoc
       to help locate it in the original file from which it was read."""

    __slots__ = ["text", "srcloc"]

    def __init__(self, text, srcloc):  # srcloc_end):
        self.text = text
        if srcloc == None:
            self.srcloc = OSrcLoc()
        else:
            self.srcloc = srcloc
        # if srcloc_end == None:
        #    self.srcloc_end = OSrcLoc()
        # else:
        #    self.srcloc_end = srcloc_end

    def __repr__(self):
        return '\"' + self.text + '\"'


class VarRef(object):
    """VarRef stores variable names, and paths, and other attribute information,
    as well as a "OSrcLoc" to keep track of the file it was defined in."""

    __slots__ = ["prefix", "descr_str", "suffix", "srcloc", "binding", "nptr"]

    def __init__(self,
                 prefix='',  # '$' or '${'
                 descr_str='',  # <- descriptor string: "cpath/category:lpath"
                 suffix='',  # '}'
                 srcloc=None,  # location in file where defined
                 binding=None,  # a pointer to a tuple storing the value
                 nptr=None):  # <- see class VarNPtr

        self.prefix = prefix  # Any text before the descriptor string goes here
        self.suffix = suffix  # Any text after the descriptor string goes here
        self.descr_str = descr_str
        if srcloc == None:  # <- Location in text file where variable appears
            self.srcloc = OSrcLoc()
        else:
            self.srcloc = srcloc

        self.binding = binding

        if nptr == None:
            self.nptr = VarNPtr()
        else:
            self.nptr = nptr

    def __lt__(self, x):
        return self.order < x.order

    # def __repr__(self):
    #    return repr((self.prefix + self.descr_str + self.suffix, srcloc))


class VarNPtr(object):
    """
    Every time a variable appears in a template, it has has a "descriptor".
    For example, consider the variable
       "$atom:CA"
    This is a string which encodes 3 pieces of information.
    1) the category name:  This is essentialy indicates the variable's type.
                           (ie "atom", in the example above)
    2) the category node:  Some TYPES have limited scope. Users can
                           specify the root node of the portion of the tree
                           in which this variable's type makes sense.
                           If this node is the root node, then that category
                           is relevant everywhere, and is not molecule or class
                           specific.  All variables have a category node, which
                           is often not explicitly defined to by the user.
                           (Category node = the root "/", in the example above.)
    3) the leaf node:      This is a node whose ".name" member matches the name
                           of a variable.  This node is created for this purpose
                           and it's position in the tree is a reflection of
                           that variable's intended scope.
                              In a molecule this "name" might be the name
                           of a type of atom, or an atom ID, or a bond type,
                           which is found in a particular molecule.
                           (Leaf node would be named "CA" in the example above.)

    The VarNPtr class is simply a 3-tuple which
    keeps these 3 pieces of data together.

    """

    __slots__ = ["cat_name", "cat_node", "leaf_node"]

    def __init__(self, cat_name='', cat_node=None, leaf_node=None):
        self.cat_name = cat_name
        self.cat_node = cat_node
        self.leaf_node = leaf_node

    # def __repr__(self):
    #    return repr((self.cat_name, self.cat_node.name, self.leaf_node.name))


class VarBinding(object):
    """ VarBinding is essentially a tuple consistng of (full_name, binding, refs):

    "self.full_name" is canonical name for this variable.  This is a string
    which specifies full path leading to the category node (beginning with '/'),
    the category name (followed by a ':'),
    as well as the leaf node (including the path leading up to it from cat_node)
    This triplet identifies the variable uniquely.

    "self.value" is the data that the variable refers to (usually a string).

    "self.refs" stores a list of VarRefs which mention the same variable
    from the various places inside various templates in the tree.

    """

    __slots__ = ["full_name", "nptr", "value", "refs", "order", "category"]

    def __init__(self,
                 full_name='',
                 nptr=None,
                 value=None,
                 refs=None,
                 order=-1,
                 category=None):
        self.full_name = full_name
        self.nptr = nptr
        self.value = value
        self.refs = refs
        self.order = order
        self.category = category

    def __lt__(self, x):
        return self.order < x.order

    def __repr__(self):
        return repr((self.full_name, self.value, self.order))


def ExtractCatName(descr_str):
    """ When applied to a VarRef's "descr_str" member,
    this function will extract the "catname" of it's corresponding
    "nptr" member.  This can be useful for error reporting.
    (I use it to insure that the user is using the correct counter
     variable types at various locations in their input files.)

    """

    ib = descr_str.find(':')
    if ib == -1:
        ib = len(descr_str)
        ia = descr_str.rfind('/')
        if ia == -1:
            ia = 0
        return descr_str[ia:ib]
    else:
        str_before_colon = descr_str[0:ib]
        ia = str_before_colon.rfind('/')
        if ia == -1:
            return str_before_colon
        else:
            return str_before_colon[ia + 1:]


def _DeleteLineFromTemplate(tmpl_list,
                            i_entry,  # index into tmpl_list
                            newline_delimiter='\n'):
    """ Delete a single line from tmpl_list.
    tmpl_list is an alternating list of VarRefs and TextBlocks.
    To identify the line, the index corresponding to one of the
    entries in the tmpl_list is used. (Usually it is a VarRef)
    The text after the preceeding newline, and the text up to the next newline
       (starting from the beginning of the current entry, if a TextBlock)
    is deleted, including any VarRef (variables) located in between.

    It returns the index corresponding to the next
    entry in the list (after deletion).

    """

    i_prev_newline = i_entry
    while i_prev_newline >= 0:
        entry = tmpl_list[i_prev_newline]
        if isinstance(entry, TextBlock):
            i_char_newline = entry.text.rfind(newline_delimiter)
            if i_char_newline != -1:  # then newline found
                # Delete the text after this newline
                entry.text = entry.text[:i_char_newline + 1]
                break
        i_prev_newline -= 1

    first_var = True
    #i_next_newline = i_entry
    i_next_newline = i_prev_newline + 1
    while i_next_newline < len(tmpl_list):
        entry = tmpl_list[i_next_newline]
        if isinstance(entry, TextBlock):
            i_char_newline = entry.text.find(newline_delimiter)
            if i_char_newline != -1:  # then newline found
                # Delete the text before this newline (including the newline)
                entry.text = entry.text[i_char_newline + 1:]
                break
        # Invoke DeleteSelf() on the first variables on this line.  This will
        # insure that it is deleted from the ttree_assignments.txt file.
        elif isinstance(entry, VarRef):
            if first_var:
                entry.nptr.leaf_node.DeleteSelf()
            first_var = False
        i_next_newline += 1

    del tmpl_list[i_prev_newline + 1: i_next_newline]
    return i_prev_newline + 1


def DeleteLinesWithBadVars(tmpl_list,
                           delete_entire_template=False,
                           newline_delimiter='\n'):
    """
    Loop through the entries in a template,
    an alternating list of TextBlocks and VarRefs (tmpl_list).
    If a VarRef points to a leaf_node which no longer exists
    (ie. no longer in the corresponding category's .bindings list).
    Then delete the line it came from from the template (tmpl_list).

    """

    out_str_list = []
    i = 0
    while i < len(tmpl_list):
        entry = tmpl_list[i]
        if isinstance(entry, VarRef):
            var_ref = entry
            var_bindings = var_ref.nptr.cat_node.categories[
                var_ref.nptr.cat_name].bindings
            # if var_ref.nptr.leaf_node not in var_bindings:
            if var_ref.nptr.leaf_node.IsDeleted():
                if delete_entire_template:
                    del tmpl_list[:]
                    return 0
                else:
                    i = _DeleteLineFromTemplate(tmpl_list,
                                                i,
                                                newline_delimiter)
            else:
                i += 1
        else:
            i += 1


def SplitTemplate(ltmpl, delim, delete_blanks=False):
    """
    Split a template "ltmpl" into a list of "tokens" (sub-templates)
    using a single delimiter string "delim".

    INPUT arguments:
    "ltmpl" should be an list of TextBlocks and VarRefs.
    "delim" should be a simple string (type str)
    "delete_blanks" should be a boolean True/False value.
                    When true, successive occurrences of the delimiter
                    should not create blank entries in the output list.

    OUTPUT:
    A list of tokens.
    Each "token" is either a TextBlock, a VarRef,
    or a (flat, 1-dimensional) list containing more than one of these objects.
    The number of "tokens" returned equals the number of times the delimiter
    is encountered in any of the TextBlocks in the "ltmpl" argument, plus one.
    (... Unless "delete_blanks" is set to True.
     Again, in that case, empty entries in this list are deleted.)

    """
    assert(type(delim) is str)
    if not hasattr(ltmpl, '__len__'):
        ltmpl = [ltmpl]

    tokens_lltmpl = []
    token_ltmpl = []
    i = 0
    while i < len(ltmpl):

        entry = ltmpl[i]
        #sys.stderr.write('ltmpl['+str(i)+'] = '+str(entry)+'\n')

        if isinstance(entry, TextBlock):
            # if hasattr(entry, 'text'):
            prev_src_loc = entry.srcloc

            tokens_str = entry.text.split(delim)

            lineno = entry.srcloc.lineno

            j = 0
            while j < len(tokens_str):
                token_str = tokens_str[j]

                delim_found = False
                if (j < len(tokens_str) - 1):
                    delim_found = True

                if token_str == '':
                    if delete_blanks:
                        if delim == '\n':
                            lineno += 1
                        if len(token_ltmpl) > 0:
                            if len(token_ltmpl) == 1:
                                tokens_lltmpl.append(token_ltmpl[0])
                            else:
                                tokens_lltmpl.append(token_ltmpl)
                        del token_ltmpl
                        token_ltmpl = []
                        j += 1
                        continue

                new_src_loc = OSrcLoc(prev_src_loc.infile, lineno)
                new_src_loc.order = prev_src_loc.order

                for c in token_str:
                    # Reminder to self:  c != delim  (so c!='\n' if delim='\n')
                    # (We keep track of '\n' characters in delimiters above.)
                    if c == '\n':
                        lineno += 1

                new_src_loc.lineno = lineno

                text_block = TextBlock(token_str,
                                       new_src_loc)

                prev_src_loc = new_src_loc

                if len(token_ltmpl) == 0:
                    if delim_found:
                        tokens_lltmpl.append(text_block)
                        del token_ltmpl
                        token_ltmpl = []
                    else:
                        token_ltmpl.append(text_block)
                else:
                    if delim_found:
                        if len(token_str) > 0:
                            token_ltmpl.append(text_block)
                            tokens_lltmpl.append(token_ltmpl)
                            del token_ltmpl
                            token_ltmpl = []
                        else:
                            assert(not delete_blanks)
                            if (isinstance(token_ltmpl[-1], VarRef)
                                and
                                ((j > 0)
                                 or
                                 ((j == len(tokens_str) - 1) and
                                  (i == len(ltmpl) - 1))
                                 )):
                                # In that case, this empty token_str corresponds
                                # to a delimiter which was located immediately
                                # after the variable name,
                                # AND
                                #   -there is more text to follow,
                                #   OR
                                #   -we are at the end of the template.
                                token_ltmpl.append(text_block)
                            if len(token_ltmpl) == 1:
                                tokens_lltmpl.append(token_ltmpl[0])
                            else:
                                tokens_lltmpl.append(token_ltmpl)
                            del token_ltmpl
                            token_ltmpl = []
                    else:
                        token_ltmpl.append(text_block)

                if (delim_found and (delim == '\n')):
                    lineno += 1

                j += 1

        elif isinstance(entry, VarRef):
            # elif hasattr(entry, 'descr_str'):
            lineno = entry.srcloc.lineno
            if ((len(token_ltmpl) == 1) and
                    isinstance(token_ltmpl[0], TextBlock) and
                    (len(token_ltmpl[0].text) == 0)):
                # special case: if the previous entry was "", then it means
                # the delimeter appeared at the end of the previous text block
                # leading up to this variable.  It separates the variable from
                # the previous text block.  It is not a text block of length 0.
                token_ltmpl[0] = entry
            else:
                token_ltmpl.append(entry)
        elif entry == None:
            token_ltmpl.append(entry)
        else:
            assert(False)

        i += 1

    # Append left over remains of the last token
    if len(token_ltmpl) == 1:
        tokens_lltmpl.append(token_ltmpl[0])
    elif len(token_ltmpl) > 1:
        tokens_lltmpl.append(token_ltmpl)
    del token_ltmpl

    return tokens_lltmpl


def SplitTemplateMulti(ltmpl, delims, delete_blanks=False):
    """
    Split a template "ltmpl" into a list of templates using a
    single one or more delimiter strings "delim_list".
    If multiple delimiter strings are provided, splitting
    begins using the first delimiter string in the list.
    Then each token in the resulting list of templates
    is split using the next delimiter string
    and so on until we run out of delimiter strings.

    "ltmpl" should be an list of TextBlocks and VarRefs.
    "delims" should be a simple string (type str) or a list of strings
    "delete_blanks" is either True or False
                    If True, then any blank entries in the resulting list of
                    tokens (sub-templates) will be deleted.

    """

    if hasattr(delims, '__len__'):  # then it hopefully is a list of strings
        delim_list = delims
    else:
        delim_list = [delims]     # then it hopefully is a string

    tokens = [ltmpl]
    for delim in delim_list:
        assert(type(delim) is str)
        tokens_il = []
        for t in tokens:
            sub_tokens = SplitTemplate(t, delim, delete_blanks)
            for st in sub_tokens:
                if hasattr(st, '__len__'):
                    if (len(st) > 0) or (not delete_blanks):
                        tokens_il.append(st)
                else:
                    tokens_il.append(st)
        tokens = tokens_il
        del tokens_il

    return tokens


def _TableFromTemplate(d, ltmpl, delimiters, delete_blanks):
    """
    See the docstring for the TableFromTemplate() function for an explanation.
    (This _TableFromTemplate() and SplitTemplate() are the workhorse functions
     for TableFromTemplate().)

    """

    output = SplitTemplateMulti(ltmpl, delimiters[d], delete_blanks[d])

    if d > 0:
        i = 0
        while i < len(output):
            output[i] = _TableFromTemplate(d - 1,
                                           output[i],
                                           delimiters,
                                           delete_blanks)
            # Delete empty LISTS?
            if (delete_blanks[d] and
                    hasattr(output[i], '__len__') and
                    (len(output[i]) == 0)):
                del output[i]
            else:
                i += 1

    return output


def TableFromTemplate(ltmpl, delimiters, delete_blanks=True):
    """
    This function can be used to split a template
    (a list containing TextBlocks and VarRefs) into a table
    into a multidimensional table, with an arbitrary number of dimensions.

    Arguments:

    ltmpl

    An alternating list of TextBlocks and VarRefs containing
    the contents of this text template.

    delimiters

    The user must supply a list or tuple of delimiters: one delimiter for
    each dimension in the table, with low-priority delimiters
    (such as spaces ' ') appearing first, and higher-priority delimiters
    (sich as newlines '\n') appearing later on in the list.
    This function will divide the entire "ltmpl" into an n-dimensional
    table.  Initially the text is split into a list of text using the
    highest-priority delimiter.  Then each entry in the resulting list is
    split into another list according to the next highest-priority delimiter.
    This continues until all of the delimiters are used up and an
    n-dimensional list-of-lists is remaining.

    delete_blanks

    The optional "delete_blanks" argument can be used to indicate whether
    or not to delete blank entries in the table (which occur as a result
    of placing two delimiters next to each other).  It should be either
    None (default), or it should be an array of booleans matching the
    size of the "delimiters" argument.  This allows the caller to customize
    the merge settings separately for each dimension (for example: to allow
    merging of whitespace within a line, without ignoring blank lines).


     ---- Details: ----

    1) Multi-character delimiters ARE allowed (like '\n\n').

    2) If a delimiter in the "delimiters" argument is not a string
    but is a tuple (or a list) of strings, then the text is split according
    to any of the delimiters in that tuple/list (starting from the last entry).
    This way, users can use this feature to split text according to multiple
    different kinds of whitespace characters (such as ' ' and '\t'), for
    example, buy setting delimiters[0] = (' ','\t').   If, additionally,
    delete_blanks[0] == True, then this will cause this function to
    divide text in without regard to whitespace on a given line (for example).

    Detailed example:

    table2D = TableFromTmplList(ltmpl,
                                delimiters = ((' ','\t'), '\n'),
                                delete_blanks = (True, False))

    This divides text in a similar way that the "awk" program does by default,
    ie, by ignoring various kinds of whitespace between text fields, but NOT
    ignoring blank lines.

    3) Any text contained in variable-names is ignored.

    """

    # Make a copy of ltmpl
    # (The workhorse function "_TableFromTemplate()" makes in-place changes to
    #  its "ltmpl" argument.  I don't want to modify "ltmpl", so I make a copy
    #  of it before I invoke "_TableFromTemplate()" on it.)

    output = [ltmpl[i] for i in range(0, len(ltmpl))]

    d = len(delimiters) - 1
    output = _TableFromTemplate(d, output, delimiters, delete_blanks)
    return output


class TemplateLexer(TtreeShlex):
    """ This class extends the standard python lexing module, shlex, adding a
    new member function (ReadTemplate()), which can read in a block of raw text,
    (halting at an (non-escaped) terminal character), and split the text into
    alternating blocks of text and variables.  (As far as this lexer is
    concerned, "variables" are simply tokens preceeded by $ or @ characters,
    and surrounded by optional curly-brackets {}.)

    """

    def __init__(self,
                 instream=None,
                 infile=None,
                 posix=False):
        TtreeShlex.__init__(self, instream, infile, posix)
        self.var_delim = '$@'  # characters which can begin a variable name
        self.var_open_paren = '{'  # optional parenthesis surround a variable
        self.var_close_paren = '}'  # optional parenthesis surround a variable
        self.newline = '\n'
        self.comment_skip_var = '#'

        #   Which characters belong in words?
        #
        # We want to allow these characters:
        #     ./$@&%^!*~`-_:;?<>[]()
        # to appear inside the tokens that TtreeShlex.get_token()
        # retrieves (TtreeShlex.get_token() is used to read class
        # names, and instance names, and variable names)
        #
        # settings.lex.wordchars+='./$@&%^!*~`-_+:;?<>[]' #Allow these chars
        #
        # Ommisions:
        # Note: I left out quotes, whitespace, comment chars ('#'), and escape
        #       characters ('\\') because they are also dealt with separately.
        #       Those characters should not overlap with settings.lex.wordchars.
        #
        # Enabling unicode support requires that we override this choice
        # by specifying "lex.wordterminators" instead of "wordchars".
        #
        # lex.wordterminators should be the (printable) set inverse of lex.wordchars
        # I'm not sure which ascii characters are NOT included in the string above
        # (We need to figure that out, and put them in settings.lex.wordterminators)
        # To figure that out, uncomment the 8 lines below:
        #
        # self.wordterminators=''
        # for i in range(0,256):
        #    c = chr(i)
        #    if c not in self.wordchars:
        #        self.wordterminators += c
        #sys.stderr.write('-------- wordterminators = --------\n')
        # sys.stderr.write(self.wordterminators+'\n')
        # sys.stderr.write('-----------------------------------\n')
        #
        # Here is the result:
        self.wordterminators = '(){|}' + \
            self.whitespace + \
            self.quotes + \
            self.operators + \
            self.escape + \
            self.commenters

        #  Note:
        # self.whitespace = ' \t\r\f\n'
        # self.quotes     = '\'"'
        # self.escape     = '\\'
        # self.commenters = '#'
        #  Note: I do not terminate on these characters: +-=*'"`
        # because they appear in the names of atom types in many force-fields.
        # Also * characters are needed for variables containing wildcards
        # in the name (which will be dealt with later).

        self.source_triggers = set(['include', 'import'])
        self.source_triggers_x = set(['import'])

    def GetSrcLoc(self):
        return OSrcLoc(self.infile, self.lineno)

    def ReadTemplate(self,
                     simplify_output=False,
                     terminators='}',
                     remove_esc_preceeding='{\\',  #explained below
                     var_terminators='{}(),', #(var_delim, spaces also included)
                     keep_terminal_char=True):
        """
           ReadTemplate() reads a block of text (between terminators)
        and divides it into variables (tokens following a '$' or '@' character)
        and raw text.  This is similar to pythons string.Template(),
        however it reads from streams (files), not strings, and it allows use
        of more complicated variable names with multiple variable delimiters
        (eg '$' and '@').
        This readline()-like member function terminates when reaching a
        user-specified terminator character character (second argument),
        or when variable (eg: "$var"$ is encountered).  The result is
        a list of variable-separated text-blocks (stored in the first
        argument).   For example, the string:
        "string with $var1 and $var2 variables.}"  contains:
                "string with ",
                $var1,
                " and ",
                $var2,
                " variables.}"
        This simplifies the final process of rendering
        (substituting text into) the text blocks later on.
            Output:
        This function returns a list of (alternating) blocks of
        text, and variable names.  Each entry in the list is either:
        1) a text block:
               Raw text is copied from the source, verbatim, along with
               some additional data (filename and line numbers), to
               help retroactively identify where the text came from
               (in case a syntax error in the text is discovered later).
               In this case, the list entry is stored as a list
               The format (TextBlock) is similar to:
                  [text_string, ((filenameA,lineBegin), (filenameB,lineEnd))],
               where the tuples, (filenameA,lineBegin) and (filenameB,lineEnd)
               denote the source file(s) from which the text was read, and
               line number at the beginning and ending of the text block.
               (This information is useful for generating helpful error
               messages.  Note that the "TtreeShlex" class allows users to
               combine multiple files transparently into one stream using
               the "source" (or "sourcehook()") member.  For this reason, it
               is possible, although unlikely, that the text-block
               we are reading could span multiple different files.)
        2) a variable (for example "$var" or "${var}"):
               In this case, the list entry is stored in the "VarRef" format
               which is essentialy shown below:
                  [[var_prefix, var_nptr, var_suffix], (filename,lineno)]
               where var_prefix and var_suffix are strings containing brackets
               and other text enclosing the variable name (and may be empty).

       As an example, we consider a file named  "datafile" which
       contains the text containing 2 text blocks and 1 variable:
               "some\n text\n before ${var}. Text after\n".
       ReadTemplate() will read this and return a list with 3 entries:
             [ ['some\n text\n before', (('datafile', 1), ('datafile', 3))],
               [['${', 'var', '}'], ('datafile', 3, 3)],
               ['Text after\n', (('datafile', 3), ('datafile', 4))] ]

        Note that while parsing the text, self.lineno counter is
        incremented whenever a newline character is encountered.
        (Also: Unlike shlex.get_token(), this function does not
        delete commented text, or insert text from other files.)

            Exceptional Cases:
        Terminator characters are ignored if they are part of a variable
        reference. (For example, the '}' in "${cat:var}", is used to denote a
        bracketed variable, and does not cause ReadTemplate() to stop reading)
           OR if they are part of a two-character escape sequence
        (for example, '}' in "\}" does not cause terminate parsing).
        In that case, the text is considered normal text.  (However the
        \ character is also stripped out.  It is also stripped out if it
        preceeds any characters in "remove_esc_preceeding", which is
        the second argument.  Otherwise it is left in the text block.)
         What is the purpose of "remove_esc_preceeding"? To force ReadTemplate()
        to remove the preceeding \ when it otherwise would not.  For example,
        we want to remove \ whenever it preceeds another \ character, so we
        include it in the remove_esc_preceeding string variable. We alse include
        '{' because we want to remove \ when it preceeds the '{' character.
        That way the \ gets deleted when it preceeds either '{' or '}'.
        (The \ character is already removed before the '}' character.)
        We want consistent behavior that people expect, so that
        "\{abc\}" -> ReadTemplate() -> "{abc}"    (instead of "\{abc}").
        In retrospect, perhaps this is a confusing way to implement this.
        
        """

        #sys.stderr.write('    ReadTemplate('+terminators+') invoked at '+self.error_leader())

        # The main loop of the parser reads only one variable at time.
        # The following variables keep track of where we are in the template.
        reading_var = False  # Are we currently reading in the name of a variable?

        prev_char_delim = False  # True iff we just read a var_delim character like '$'
        # True iff we just read a (non-escaped) esc character '\'
        escaped_state = False
        # True iff we are in a region of text where vars should be ignored
        commented_state = False
        var_paren_depth = 0  # This is non-zero iff we are inside a
        # bracketed variable's name for example: "${var}"
        var_terminators += self.whitespace + self.newline + self.var_delim

        tmpl_list = []  # List of alternating tuples of text_blocks and
        # variable names (see format comment above)
        # This list will be returned to the caller.

        # sys.stderr.write('report_progress='+str(report_progress))

        prev_filename = self.infile
        prev_lineno = self.lineno
        var_prefix = ''
        var_descr_plist = []
        var_suffix = ''
        text_block_plist = []

        done_reading = False

        while not done_reading:

            terminate_text = False
            terminate_var = False
            #delete_prior_escape = False

            nextchar = self.read_char()

            #sys.stderr.write('    ReadTemplate() nextchar=\''+nextchar+'\' at '+self.error_leader()+'  esc='+str(escaped_state)+', pvar='+str(prev_char_delim)+', paren='+str(var_paren_depth))

            # Count newlines:
            if nextchar in self.newline:
                commented_state = False
                self.lineno += 1

            elif ((nextchar in self.comment_skip_var) and
                  (not escaped_state)):
                commented_state = True

            # Check for end-of-file:
            if nextchar == '':

                if escaped_state:
                    raise InputError('Error: in ' + self.error_leader() + '\n\n'
                                     'File terminated immediately following an escape character.')
                    terminate_var = True
                else:
                    terminate_text = True

                done_reading = True

            # --- Now process the character: ---

            # What we do next depends on which "mode" we are in.
            #  If we are reading a regular text block (reading_var == False),
            #   then we keep appending characters onto the end of "text_block",
            #   checking for terminal characters, or variable delimiters.
            #  If we are reading a variable name (reading_var == True),
            #   then we append characters to the end of "var_descr_plist[]",
            #   checking for variable terminator characters, as well as
            #   parenthesis (some variables are surrounded by parenthesis).

            elif reading_var:

                if nextchar in terminators:
                    #sys.stdout.write('   ReadTemplate() readmode found terminator.\n')
                    if escaped_state:
                        # In this case, the '\' char was only to prevent terminating
                        # string prematurely, so delete the '\' character.
                        #delete_prior_escape = True
                        del var_descr_plist[-1]
                        var_descr_plist.append(nextchar)
                        #escaped_state = False
                    elif not ((var_paren_depth > 0) and
                              (nextchar in self.var_close_paren)):
                        terminate_var = True
                        done_reading = True

                if nextchar in self.var_open_paren:  # eg: nextchar == '{'
                    #sys.stdout.write('   ReadTemplate() readmode found {\n')
                    if escaped_state:
                        var_descr_plist.append(nextchar)
                        #escaped_state = False
                    else:
                        # "${var}" is a valid way to refer to a variable
                        if prev_char_delim:
                            var_prefix += nextchar
                            var_paren_depth = 1
                        # "${{var}}" is also a valid way to refer to a variable,
                        # (although strange), but "$va{r}" is not.
                        # Parenthesis (in bracketed variable names) must
                        # immediately follow the '$' character (as in "${var}")
                        elif var_paren_depth > 0:
                            var_paren_depth += 1
                            var_descr_plist.append(nextchar)

                elif nextchar in self.var_close_paren:
                    #sys.stdout.write('   ReadTemplate() readmode found }.\n')
                    if escaped_state:
                        # In this case, the '\' char was only to prevent
                        # interpreting '}' as a variable suffix,
                        # delete_prior_escape=True  #so skip the '\' character
                        del var_descr_plist[-1]
                        var_descr_plist.append(nextchar)
                        #escaped_state = False
                    else:
                        if var_paren_depth > 0:
                            var_paren_depth -= 1
                            if var_paren_depth == 0:
                                var_suffix = nextchar
                                terminate_var = True
                            else:
                                var_descr_plist.append(nextchar)

                elif nextchar in var_terminators:
                    #sys.stdout.write('   ReadTemplate() readmode found var_terminator \"'+nextchar+'\"\n')
                    if (escaped_state or (var_paren_depth > 0)):
                        # In that case ignore the terminator
                        # and append it to the variable name
                        if escaped_state:
                            # In this case, the '\' char was only to prevent
                            # interpreting nextchar as a variable terminator
                            # delete_prior_escape = True # so skip the '\'
                            #                            # character
                            del var_descr_plist[-1]
                            #escaped_state = False
                        var_descr_plist.append(nextchar)
                    else:
                        terminate_var = True

                elif nextchar in self.var_delim:   # such as '$'
                    #sys.stdout.write('   ReadTemplate() readmode found var_delim.\n')
                    if escaped_state:
                        # In this case, the '\' char was only to prevent
                        # interpreting '$' as a new variable name
                        # delete_prior_escape = True # so skip the '\'
                        # character
                        del var_descr_plist[-1]
                        var_descr_plist.append(nextchar)
                        #escaped_state = False
                    else:
                        prev_var_delim = True
                        # Then we are processing a new variable name
                        terminate_var = True
                else:
                    var_descr_plist.append(nextchar)
                    prev_char_delim = False

            else:  # begin else clause for "if reading_var:"

                # Then we are reading a text_block

                if nextchar in terminators:
                    if escaped_state:
                        # In this case, the '\' char was only to prevent terminating
                        # string prematurely, so delete the '\' character.
                        #delete_prior_escape = True
                        del text_block_plist[-1]
                        text_block_plist.append(nextchar)
                    elif commented_state:
                        text_block_plist.append(nextchar)
                    else:
                        terminate_text = True
                        done_reading = True

                elif nextchar in self.var_delim:   # such as '$'
                    if escaped_state:
                        # In this case, the '\' char was only to prevent
                        # interpreting '$' as a variable prefix.
                        # delete_prior_escape=True  #so delete the '\'
                        # character
                        del text_block_plist[-1]
                        text_block_plist.append(nextchar)
                    elif commented_state:
                        text_block_plist.append(nextchar)
                    else:
                        prev_char_delim = True
                        reading_var = True
                        # NOTE TO SELF: IN THE FUTURE, USE GetVarName(self)
                        # TO PARSE TEXT ASSOCIATED WITH A VARIABLE
                        # THIS WILL SIMPLIFY THE CODE AND ENSURE CONSISTENCY.
                        var_paren_depth = 0
                        terminate_text = True
                else:
                    text_block_plist.append(nextchar)
                    # TO DO: use "list_of_chars.join()" instead of '+='
                    prev_char_delim = False  # the previous character was not '$'

            # Now deal with "remove_esc_preceeding".  (See explanation above.)
            if escaped_state and (nextchar in remove_esc_preceeding):
                if reading_var:
                    #sys.stdout.write('   ReadTemplate: var_descr_str=\''+''.join(var_descr_plist)+'\'\n')
                    assert(var_descr_plist[-2] in self.escape)
                    del var_descr_plist[-2]
                else:
                    #sys.stdout.write('   ReadTemplate: text_block=\''+''.join(text_block_plist)+'\'\n')
                    assert(text_block_plist[-2] in self.escape)
                    del text_block_plist[-2]

            if terminate_text:
                #sys.stdout.write('ReadTemplate() appending: ')
                # sys.stdout.write(text_block)

                # tmpl_list.append( [text_block,
                #                   ((prev_filename, prev_lineno),
                #                    (self.infile, self.lineno))] )

                if simplify_output:
                    tmpl_list.append(''.join(text_block_plist))
                else:
                    tmpl_list.append(TextBlock(''.join(text_block_plist),
                                               OSrcLoc(prev_filename, prev_lineno)))
                    #, OSrcLoc(self.infile, self.lineno)))
                if not done_reading:
                    # The character that ended the text block
                    # was a variable delimiter (like '$'), in which case
                    # we should put it (nextchar) in the variable's prefix.
                    var_prefix = nextchar
                else:
                    var_prefix = ''
                var_descr_plist = []
                var_suffix = ''
                prev_filename = self.infile
                prev_lineno = self.lineno
                del text_block_plist
                text_block_plist = []
                # gc.collect()

            elif terminate_var:
                # Print an error if we terminated in the middle of
                # an incomplete variable name:
                if prev_char_delim:
                    raise InputError('Error: near ' + self.error_leader() + '\n\n'
                                     'Null variable name.')
                if var_paren_depth > 0:
                    raise InputError('Error: near ' + self.error_leader() + '\n\n'
                                     'Incomplete bracketed variable name.')

                var_descr_str = ''.join(var_descr_plist)

                # Now check for variable format modifiers,
                # like python's ".rjust()" and ".ljust()".
                # If present, then put these in the variable suffix.
                if ((len(var_descr_plist) > 0) and (var_descr_plist[-1] == ')')):
                    #i = len(var_descr_plist)-1
                    # while i >= 0:
                    #    if var_descr_plist[i] == '(':
                    #        break
                    #    i -= 1
                    i = var_descr_str.rfind('(')
                    if (((i - 6) >= 0) and
                        ((var_descr_str[i - 6:i] == '.rjust') or
                         (var_descr_str[i - 6:i] == '.ljust'))):
                        var_suffix = ''.join(
                            var_descr_plist[i - 6:]) + var_suffix
                        #var_descr_plist = var_descr_plist[:i-6]
                        var_descr_str = var_descr_str[:i - 6]

                # Process any special characters in the variable name
                var_descr_str = EscCharStrToChar(var_descr_str)

                # tmpl_list.append( [[var_prefix, var_descr_str, var_suffix],
                #                   (self.infile, self.lineno)] )
                if simplify_output:
                    tmpl_list.append(var_prefix + var_descr_str + var_suffix)
                else:
                    tmpl_list.append(VarRef(var_prefix, var_descr_str, var_suffix,
                                            OSrcLoc(self.infile, self.lineno)))

                # if report_progress:
                #sys.stderr.write('  parsed variable '+var_prefix+var_descr_str+var_suffix+'\n')

                #sys.stdout.write('ReadTemplate() appending: ')
                #sys.stderr.write(var_prefix + var_descr_str + var_suffix)

                del var_descr_plist
                del var_descr_str

                prev_filename = self.infile
                prev_lineno = self.lineno
                var_prefix = ''
                var_descr_plist = []
                var_suffix = ''
                # Special case: Variable delimiters like '$'
                #               terminate the reading of variables,
                #               but they also signify that a new
                #               variable is being read.
                if nextchar in self.var_delim:
                    # Then we are processing a new variable name
                    prev_var_delim = True
                    reading_var = True
                    # NOTE TO SELF: IN THE FUTURE, USE GetVarName(self)
                    # TO PARSE TEXT ASSOCIATED WITH A VARIABLE
                    # THIS WILL SIMPLIFY THE CODE AND ENSURE CONSISTENCY.
                    var_paren_depth = 0
                    var_prefix = nextchar

                elif nextchar in self.var_close_paren:
                    del text_block_plist
                    text_block_plist = []
                    # gc.collect()
                    prev_var_delim = False
                    reading_var = False

                else:
                    # Generally, we don't want to initialize the next text block
                    # with the empty string.  Consider that whatever character
                    # caused us to stop reading the previous variable and append
                    # it to the block of text that comes after.
                    del text_block_plist
                    text_block_plist = [nextchar]
                    # gc.collect()
                    prev_var_delim = False
                    reading_var = False

            # If we reached the end of the template (and the user requests it),
            # then the terminal character can be included in the list
            # of text_blocks to be returned to the caller.
            if done_reading and keep_terminal_char:
                #sys.stdout.write('ReadTemplate() appending: \''+nextchar+'\'\n')
                # Here we create a new text block which contains only the
                # terminal character (nextchar).
                # tmpl_list.append( [nextchar,
                #                   ((self.infile, self.lineno),
                #                    (self.infile, self.lineno))] )
                if simplify_output:
                    tmpl_list.append(nextchar)
                else:
                    tmpl_list.append(TextBlock(nextchar,
                                               OSrcLoc(self.infile, self.lineno)))
                    #, OSrcLoc(self.infile, self.lineno)))

            if escaped_state:
                escaped_state = False
            else:
                if nextchar in self.escape:
                    escaped_state = True

        #sys.stderr.write("*** TMPL_LIST0  = ***", tmpl_list)
        return tmpl_list  # <- return value stored here

    def GetParenExpr(self, prepend_str='', left_paren='(', right_paren=')'):
        """ GetParenExpr() is useful for reading in strings
            with nested parenthesis and spaces.
            This function can read in the entire string:

              .trans(0, 10.0*sin(30), 10.0*cos(30))

            (Because I was too lazy to write this correctly...)
            Spaces are currently stripped out of the expression.
            (...unless surrounded by quotes) The string above becomes:

              ".trans(0,10.0*sin(30),10.0*cos(30))"

            Sometimes the caller wants to prepend some text to the beginning
            of the expression (which may contain parenthesis).  For this
            reason, an optional first argument ("prepend_str") can be
            provided.  By default it is empty.

        """

        src_loc_begin = SrcLoc(self.infile, self.lineno)
        orig_wordterm = self.wordterminators
        self.wordterminators = self.wordterminators.replace(
            left_paren, '').replace(right_paren, '')

        token = self.get_token()
        if ((token == '') or
                (token == self.eof)):
            return prepend_str

        expr_str = prepend_str + token

        # if (expr_str.find(left_paren) == -1):
        #    raise InputError('Error near or before '+self.error_leader()+'\n'
        #                     'Expected an open-paren (\"'+prepend_str+left_paren+'\") before this point.\n')
        #    return expr_str

        paren_depth = expr_str.count(left_paren) - expr_str.count(right_paren)
        while ((len(expr_str) == 0) or (paren_depth > 0)):
            token = self.get_token()
            if ((type(token) is not str) or
                    (token == '')):
                raise InputError('Error somewhere between ' +
                                 self.error_leader(src_loc_begin.infile,
                                                   src_loc_begin.lineno)
                                 + 'and ' + self.error_leader() + '\n'
                                 'Invalid expression: \"' + expr_str[0:760] + '\"')
            expr_str += token
            paren_depth = expr_str.count(
                left_paren) - expr_str.count(right_paren)
        if (paren_depth != 0):
            raise InputError('Error somewhere between ' +
                             self.error_leader(src_loc_begin.infile,
                                               src_loc_begin.lineno)
                             + 'and ' + self.error_leader() + '\n'
                             'Invalid expression: \"' + expr_str[0:760] + '\"')
        self.wordterminators = orig_wordterm
        return expr_str


if __name__ == '__main__':
    if len(sys.argv) == 1:
        lexer = TtreeShlex()
    else:
        file = sys.argv[1]
        lexer = TtreeShlex(open(file), file)
    while 1:
        tt = lexer.get_token()
        if tt:
            sys.stderr.write("Token: " + repr(tt))
        else:
            break
