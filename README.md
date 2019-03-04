[![Build Status](https://travis-ci.org/jewettaij/moltemplate.svg?branch=master)](https://travis-ci.org/jewettaij/moltemplate.svg?branch=master)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Moltemplate
===========

##  Description

Moltemplate is a *general* cross-platform text-based molecule builder for **LAMMPS** and **ESPResSo**.  Moltemplate was intended for building custom coarse-grained molecular models, but it can be used to prepare realistic all-atom simulations as well.  It currently supports the **OPLS**, **COMPASS**, **AMBER**(GAFF,GAFF2), **MARTINI**, **SDK**, **LOPLS**(2015), **EFF**, and **TraPPE**(1998) force fields, numerous solvent models, and includes approximately 40 examples.  (New force fields and examples are added continually by users.)

## Typical usage

    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt

## Web page

Documentation, examples, and supporting code can be downloaded at:

http://www.moltemplate.org

## Requirements

Moltemplate requires the Bourne-shell, and a recent version of python
(2.7, 3.0 or higher), and can run on OS X, linux, or windows. (...if a
suitable shell environment has been installed.  See below.)


## INSTALLATION INSTRUCTIONS

This directory should contain 3 folders:

    moltemplate/                  <-- source code and force fields
    doc/                          <-- the moltemplate reference manual
    examples/                     <-- examples built with moltemplate

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
    python3 -m venv venv     #(or "virtualenv venv" if using python2)
    source venv/bin/activate
    pip install .
    #(now do something useful with moltemplate...)

(You will have to "run source ~/moltemplate/venv/bin/activate" beforehand every time you want to run moltemplate.
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

If you use the bash shell, typically you would edit your
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files
to contain the following lines:

    export PATH="$PATH:$HOME/moltemplate/moltemplate"
    export PATH="$PATH:$HOME/moltemplate/moltemplate/scripts"

If you use the tcsh shell, typically you would edit your
`~/.login`, `~/.cshrc`, or `~/.tcshrc` files to contain the following lines:

    setenv PATH "$PATH:$HOME/moltemplate/moltemplate"
    setenv PATH "$PATH:$HOME/moltemplate/moltemplate/scripts"

After making these changes, you may need to start a new terminal (shell) for the changes to take effect.  If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)


### WINDOWS installation suggestions

You can install both moltemplate and LAMMPS in windows, but you will first need to install the BASH shell environment on your computer.  I recommend installing [virtualbox](https://www.virtualbox.org) in windows together with a (debian-based) linux distribution with a lightweight desktop such as [xubuntu](https://xubuntu.org).  Alternatively, if you are using Windows 10 or later, you can try installing the "Windows Subsystem for Linux (WSL)", as explained
[here](https://solarianprogrammer.com/2017/04/15/install-wsl-windows-subsystem-for-linux/)
and
[here](https://msdn.microsoft.com/en-us/commandline/wsl/faq),
or
[Hyper-V](https://blogs.windows.com/buildingapps/2018/09/17/run-ubuntu-virtual-machines-made-even-easier-with-hyper-v-quick-create/).
Otherwise, if you are using an older version of windows, try installing
[CYGWIN](https://www.cygwin.com/) instead.

To use LAMMPS and moltemplate, you will also need to install (and learn how to use) a text editor.  (Word, Wordpad, and Notepad will not work.)  Popular free text editors which you can safely install and run from within the WSL terminal include: **nano**, **ne**, **emacs**, **vim**, and **jove**.
If you are **NOT using WSL**, then you can also use a graphical text editor which might be easier for a beginner to use.  Such editors include Atom, Sublime, Notepad++, VSCode, and the graphical versione of emacs and vim.  ([Don't use these editors to edit files within the WSL environment.](https://www.reddit.com/r/bashonubuntuonwindows/comments/6bu1d1/since_we_shouldnt_edit_files_stored_in_wsl_with/))

## License

Moltemplate is available under the terms of the open-source 3-clause BSD
license.  (See "[LICENSE.md](./LICENSE.md)")
