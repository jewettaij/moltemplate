Moltemplate Installation
====================

Before you install moltemplate,
[you must install BASH and python](#OS-specific-suggestions).

Once you have met those requirements,
there are **two ways** to install moltemplate:
1) using pip
2) editing your .bashrc file

Use one method or the other *(not both)*.


## Installation method 1: *Using pip*

If you are familiar with pip, then run the following command from
within the directory where this README file is located:
```
pip3 install .              # (or "pip", if that fails)
```
If that fails (with both "pip" and "pip3"), then try this instead:
```
pip3 install . --user       # (or "pip", if that fails)
```
This will install moltemplate for a single user.
If you are on a shared computer and you want to install moltemplate
system-wide, then use:
```
sudo pip3 install .         # (or "pip", if that fails)
```
Later, you can uninstall moltemplate using:
```
pip3 uninstall moltemplate
# (use "pip" and/or prepend "sudo" if you did that earlier)
```


### Troubleshooting (method 1)

Again, if you get the error "command not found",
try using "pip" instead of "pip3".


### Details (method 1)

When users have difficulty installing moltemplate, it is usually
for one of these reasons:
 - They have not installed BASH (a.k.a. Xcode, WSL, linux, etc...).
 - The computer has an old version of python (before python 3).
 - pip and pip3 are missing. (This is very common.)
 - There are multiple conflicting versions of python and pip installed.

Once these issues are fixed, installing moltemplate is easy.

Unfortunately, there are several different versions of "python" and "pip",
and they go by different names whose meaning is not consistent.

Furthermore, many computers come with python installed
but lack "pip" (or "pip3").  In that case, you can use the official
package management system on your computer to install "pip3" (or "pip").
Alternatively, you can try installing the
[anaconda version of python](https://anaconda.com).
*Anaconda python* is recommended because it automatically
includes *numpy*, and it will also automatically update your
[PATH](http://www.linfo.org/path_env_var.html)
for you so that *pip* works correctly.
*(Note that if you installed anaconda then use "pip",
not "pip3" in the commands above.)*

Some computers have both "pip" and "pip3" installed.
In that case, "pip" and "pip3" may refer to different programs.
This usually happens when there are multiple conflicting versions of python
installed on your computer, and it is the most frequent reason for moltemplate
installation failure.
You must choose the version ("pip" or "pip3") which corresponds to the
version of python that moltemplate is using.  To do that, enter:
```
which python3
which pip3
```
If this returns "/usr/bin/python3" and "/usr/bin/pip3"
then it is probably safe to use "pip3" in the commands above
(because in this example, both "python3" and "pip3" share the same path).

If "which python3" returns something ending in "anaconda3/bin/python3"
then try using "pip" (instead of "pip3") in the commands above.

On the other hand, if "which python3" fails to return anything, try entering:
```
which python
which pip
```
If this returns "/usr/bin/python" and "/usr/bin/pip"
then it is probably safe to use "pip" in the commands above.
In that case, try entering:
```
python --version
```
If this returns "Python 2.7.13" or something lower, I recommend upgrading
to a newer version of python.  (Python versions below 3 are no longer supported
by "pip" and may cease to work in the future.
Anaconda comes with the latest stable version of python
and you can install it with or without admin privileges.)


#### *Optional: Use a python virtual environment*

Once you have *python* and *pip* (or *python3* and *pip3*) installed,
it's never a bad idea to install moltemplate into a temporary
python "virtual environment".
In a virtual environment, it should not be necessary to use "sudo" or "--user"
to get around permissions issues that sometimes occur when using *pip*.
(You can also uninstall moltemplate cleanly simply by deleting
the directory that stores the virtual environment where it was installed.
That directory is named "venv" in the example below.)
Although using a virtual python environment should not be necessary,
if you are curious how to do it then try downloading moltemplate
to "~/moltemplate", and run these commands:

```
cd ~/moltemplate
python3 -m venv venv     #(or "virtualenv venv" if that fails)
# This will create a local directory named "venv".
# You must "activate" your environment before use:
source venv/bin/activate
# Then install moltemplate into this environment:
pip3 install .           #(or "pip" if that fails)
# (Now do something useful with moltemplate...)
```

Note that if you use a virtual environment, you will have to enter
"source ~/moltemplate/venv/bin/activate" into a terminal
before you use that terminal to run moltemplate.
Virtual environments are
[explained here](https://docs.python.org/3/tutorial/venv.html)


If all this fails, then try installing moltemplate by manually updating your
\$PATH environment variable.  Instructions for doing that are included below.
(Either way, you must have a working version of python or python3 installed.)



## Installation method 2: *Editing .bashrc*

Alternatively, you can edit your $PATH environment variable manually to
include the subdirectory where the "moltemplate.sh" script is located,
as well as the subdirectory where most of the python scripts are located.
Suppose the directory with this README file is named "moltemplate"
and is located in your home directory:

If you use the *BASH* shell, typically you would edit your
`~/.bashrc` file and add the following lines to the end of the file:

```
export PATH="$PATH:$HOME/moltemplate/moltemplate"
export PATH="$PATH:$HOME/moltemplate/moltemplate/scripts"
```

(Some people prefer to put these lines in their `~/.profile`,
 or `~/.bash_profile` files instead `~/.bashrc`.  This should work also.)
If you use the *TCSH* shell, typically you would edit your
`~/.cshrc`, `~/.tcshrc` or `~/.login` files to contain the following lines:

```
setenv PATH "$PATH:$HOME/moltemplate/moltemplate"
setenv PATH "$PATH:$HOME/moltemplate/moltemplate/scripts"
```

After making these changes, you may need to start a new terminal (shell) for the changes to take effect.  If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)

*(Warning:
Do not install moltemplate this way if you are using "vipster",
"cellpack2moltemplate", or other software that has a moltemplate python
dependency.  In order to be able to be able to run "import moltemplate"
within python, as those programs do, moltemplate must be installed using
pip or pip3.)*



# OS specific suggestions

## Linux installation suggestions

Most popular linux distributions include all of the needed prerequisites,
except for "pip3" (or "pip").  Enter "pip3" into the terminal.  If
it returns "command not found" then try again with "pip".
If that also fails, then use the package management system for your linux
distribution (eg "apt", "yum", ...) to install "pip3" or "pip".
(Eg. for ubuntu linux, use "sudo apt install python3-pip".)

Alternatively, you can install the
[anaconda version of python](https://anaconda.com)
in linux, windows, and macOS.
Anaconda python includes pip and numpy by default.
(Numpy is strongly recommended, but is not required.)


## Mac installation suggestions

Although I have never tried this, the BASH shell, python, and pip
prerequisites can be installed using Homebrew (or XCODE).
However the anaconda version of python and pip is also available for mac OS,
and it may be more reliable than the python version
included with XCODE or Homebrew.


## WINDOWS installation suggestions

You can install both moltemplate and LAMMPS in windows, but
you will first need to install the BASH shell environment on
your computer.  I recommend installing either
[Windows Subsystem for Linux (WSL2)](https://docs.microsoft.com/en-us/windows/wsl/install-win10),
***or***
[virtualbox](https://www.virtualbox.org)
(In the later case, you will also need to install a linux distribution,
preferably with a lightweight
desktop such as [xubuntu](https://xubuntu.org).)
Alternatively, you can try
[Hyper-V](https://www.nakivo.com/blog/run-linux-hyper-v/)
or (if you have an older version of windows)
[CYGWIN](https://www.cygwin.com/).

To use LAMMPS and moltemplate, you will also need to install (and learn how to
use) a (unix-style) text editor.  (Word, Wordpad, and Notepad will not work.)
Popular graphical text editors
include **Atom**, **Sublime**, **Notepad++**, and **VSCode**.
Older, non-graphical programs include **vim**, **emacs**,
**nano**, **ne**, and **jove**.

WSL and virtualbox are virtual machines that allow you to run an
alternate operating system from within windows.
In this case that operating system is linux.  The linux operating system
includes the BASH shell and python (which moltemplate needs), as well as
compilers such as g++ (which can be used for compiling LAMMPS from source),
as well as the gdb debugger (which can be useful for
understanding why LAMMPS is behaving strangely).
It creates an ideal environment for running LAMMPS and moltemplate.
WSL and virtualbox also create an alternate filesystem inside windows where
the linux operating system is stored.  Software (like moltemplate and LAMMPS)
that you install there can access the files in that filesystem.
If you **are using WSL or WSL2**, then you should
[use caution when using windows programs to edit your linux files](https://devblogs.microsoft.com/commandline/do-not-change-linux-files-using-windows-apps-and-tools/).
Consequently, it might be safer to restrict yourself to using text editors
which you have installed and can run from within the WSL or WSL2 environment.
