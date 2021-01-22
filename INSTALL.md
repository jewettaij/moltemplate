Moltemplate Installation Instructions
===============

There are two ways to install moltemplate:

## Installation using pip

If you are familiar with pip, then run the following command from within the directory where this README file is located:

    pip install . --user

This will install moltemplate for a single user.  On a shared computer, to install moltemplate system-wide, use:

    sudo pip install .

Make sure that your default pip install bin directory is in your PATH.  (This is usually something like ~/.local/bin/ or ~/anaconda3/bin/.  If you have installed anaconda, this will be done for you automatically.)  Later, you can uninstall moltemplate using:

    pip uninstall moltemplate

If you continue to run into difficulty, try installing moltemplate into a temporary virtual environment by installing "*virtualenv*", downloading moltemplate (to "~/moltemplate" in the example below), and running these commands:

    cd ~/moltemplate
    python -m venv venv     #(or "virtualenv venv" if using python2)
    source venv/bin/activate
    pip install .
    #(now do something useful with moltemplate...)

(You will have to enter "source ~/moltemplate/venv/bin/activate"
 into a terminal beforehand every time you want to run moltemplate.
Virtual environments are
[explained here](https://docs.python.org/3/tutorial/venv.html)
If all this fails, then try installing moltemplate by manually updating your
\$PATH environment variable.  Instructions for doing that are included below.


## Manual installation:

Alternatively, you can edit your $PATH environment variable manually to 
include the subdirectory where the "moltemplate.sh" script is located,
as well as the subdirectory where most of the python scripts are located.
Suppose the directory with this README file is named "moltemplate"
and is located in your home directory:

If you use the *BASH* shell, typically you would edit your
`~/.bashrc` file (or `~/.profile`, `~/.bash_profile` files)
to contain the following lines:

    export PATH="$PATH:$HOME/moltemplate/moltemplate"
    export PATH="$PATH:$HOME/moltemplate/moltemplate/scripts"

If you use the *TCSH* shell, typically you would edit your
`~/.cshrc`, `~/.tcshrc` or `~/.login` files to contain the following lines:

    setenv PATH "$PATH:$HOME/moltemplate/moltemplate"
    setenv PATH "$PATH:$HOME/moltemplate/moltemplate/scripts"

After making these changes, you may need to start a new terminal (shell) for the changes to take effect.  If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)

*(Warning:
Do not install moltemplate this way if you are using "vipster",
"cellpack2moltemplate", or other software that has a moltemplate python
dependency.  In order to be able to be able to run "import moltemplate"
within python, as these programs do, moltemplate must be installed using
pip or setuptools.)*


### WINDOWS installation suggestions

You can install both moltemplate and LAMMPS in windows, but you will first need to install the BASH shell environment on your computer.  I recommend installing [virtualbox](https://www.virtualbox.org) in windows together with a (debian-based) linux distribution with a lightweight desktop such as [xubuntu](https://xubuntu.org).  Alternatively, if you are using Windows 10 or later, you can try installing the
[Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl)
*(which is text only)*
or
[Hyper-V](https://www.nakivo.com/blog/run-linux-hyper-v/).
Otherwise, if you are using an older version of windows, try installing
[CYGWIN](https://www.cygwin.com/) instead.

To use LAMMPS and moltemplate, you will also need to install (and learn
how to use) a text editor.  (Word, Wordpad, and Notepad will not work.)
If you are **NOT using WSL**, then you can use popular graphical text editors
such as Atom, Sublime, Notepad++, VSCode,
and the graphical versione of emacs and vim.
([Don't use these editors to edit files within the WSL environment.](https://www.reddit.com/r/bashonubuntuonwindows/comments/6bu1d1/since_we_shouldnt_edit_files_stored_in_wsl_with/))
If you **ARE using WSL** then you are restricted to using non-graphical text
editors which you can safely install and run from within the WSL terminal.
These include: **nano**, **ne**, **emacs** (the text version),
**vim** (the text version), and **jove**.
