Moltemplate
===========

##  Description

Moltemplate is a cross-platform text-based molecule builder for LAMMPS. 

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

    pip install .

Make sure that your default pip install bin directory is in your PATH.  (This is usually something like ~/.local/bin/ or ~/anaconda3/bin/.  If you have installed anaconda, this will be done for you automatically.)  Later, you can uninstall moltemplate using:

    pip uninstall moltemplate

Instructions for editing your PATH are included below.  

## Manual Installation method:

Alternatively, you can edit your PATH variable manually to include
the subdirectory where the moltemplate.sh script is located.
Suppose the directory with this README file is named ``moltemplate''
and is located in your home directory:

If you use the bash shell, typically you would edit your 
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files 
to contain the following line:
    export PATH="$PATH:$HOME/moltemplate/moltemplate/scripts"
If you use the tcsh shell, typically you would edit your 
`~/.login`, `~/.cshrc`, or `~/.tcshrc` files to contain the following lines:
    setenv PATH "$PATH:$HOME/moltemplate/moltemplate/scripts"
After making these changes, you may need to start a new terminal (shell) for the changes to take effect.  If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)


### WINDOWS installation suggestions

You can install both moltemplate and LAMMPS in windows, but you will first need to install the BASH shell environment on your computer.  If you are using Windows 10 or later, try installing "Bash on Ubuntu on Windows":

https://msdn.microsoft.com/en-us/commandline/wsl/install_guide

   If you are using an older version of windows, try following the tutorial written by Yanqing Fu:

https://sourceforge.net/p/lammps/mailman/message/32599824/

To use LAMMPS and moltemplate, you will also need to install (and learn how to use) a text editor.  (Word, Wordpad, and Notepad will probably not work.)  Popular text editors include: nano, emacs, vim, gedit, xemacs, vi, jove.  A longer list can be found here:
https://en.wikipedia.org/wiki/List_of_text_editors

## License

Moltemplate is available under the terms of the open-source 3-clause BSD 
license.  (See `LICENSE.txt`.)

